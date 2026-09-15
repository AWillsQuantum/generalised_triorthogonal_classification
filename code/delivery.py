"""Verify the two deliveries or reconstruct the byte-identical combined release.

Only Python's standard library is required. ZIP extraction leaves inner Zstandard
space and assignment streams compressed. No network access is performed.
"""

import argparse
from contextlib import ExitStack
import hashlib
import json
from pathlib import Path, PurePosixPath
import shutil
import stat
import zipfile

from verify_release_integrity import verify as verify_manifest


LAYOUT = "delivery_layout.json"
ORIGINAL = "provenance/ORIGINAL_MANIFEST.json"
SCHEMA = "triorthogonal-two-delivery-layout-v2"


def safe_name(name):
    if not isinstance(name, str):
        raise ValueError("A path must be a string")
    path = PurePosixPath(name)
    if (not name or name == "." or path.is_absolute() or str(path) != name
            or ".." in path.parts or ":" in name or "\\" in name or "\x00" in name):
        raise ValueError("Unsafe relative path: " + repr(name))
    return name


def local_file(root, name):
    safe_name(name)
    path = root / name
    for parent in (path, *path.parents):
        if parent == root:
            break
        if parent.is_symlink():
            raise ValueError("Symbolic links are not permitted: " + name)
    if not path.resolve().is_relative_to(root.resolve()):
        raise ValueError("Path leaves delivery: " + name)
    return path


def stream_digest(stream):
    digest = hashlib.sha256()
    size = 0
    while block := stream.read(1024 * 1024):
        digest.update(block)
        size += len(block)
    return size, digest.hexdigest()


def catalogue_summary(raw):
    data = json.loads(raw)
    if (data.get("maximum_protocol_length") != 54 or data.get("minimum_distance") != 3
            or data.get("equivalence") != "CNOT+S"
            or data.get("intrinsic_outputs_only") is not True
            or data.get("pareto_partition") != ["output_id", "exact d_Z"]
            or data.get("pareto_objectives") != ["n", "S"]
            or len(data["protocols"]) != 74 or len(data["outputs"]) != 62):
        raise ValueError("Unexpected frontier scope or counts")
    if not all(row.get("representative_gate") for row in data["outputs"]):
        raise ValueError("A representative gate is missing")
    return {"protocols": 74, "outputs": 62, "sha256": hashlib.sha256(raw).hexdigest()}


class Delivery:
    def __init__(self, github, data):
        self.roots = {"github": Path(github).resolve(), "data": Path(data).resolve()}
        raw = local_file(self.roots["github"], LAYOUT).read_bytes()
        if raw != local_file(self.roots["data"], LAYOUT).read_bytes():
            raise ValueError("Delivery layouts differ")
        self.layout = json.loads(raw)
        if self.layout.get("schema") != SCHEMA:
            raise ValueError("Unknown delivery layout")
        original = local_file(self.roots["github"], ORIGINAL).read_bytes()
        if hashlib.sha256(original).hexdigest() != self.layout["original_manifest_sha256"]:
            raise ValueError("Original manifest does not match the layout")
        if original != local_file(self.roots["data"], "ORIGINAL_MANIFEST.json").read_bytes():
            raise ValueError("Original manifests differ")
        manifest = json.loads(original)
        if (manifest.get("schema") != "triorthogonal-release-file-manifest-v1"
                or manifest.get("algorithm") != "sha256"
                or manifest.get("excluded_files") != ["MANIFEST.json"]):
            raise ValueError("Unknown original manifest")
        self.records = {}
        for record in manifest["files"]:
            name = safe_name(record["path"])
            if name in self.records or name == "MANIFEST.json":
                raise ValueError("Repeated original path")
            if (type(record["bytes"]) is not int or record["bytes"] < 0
                    or len(record["sha256"]) != 64
                    or any(c not in "0123456789abcdef" for c in record["sha256"])):
                raise ValueError("Invalid original file metadata")
            self.records[name] = record
        if (manifest["file_count"] != len(self.records)
                or manifest["payload_bytes"] != sum(r["bytes"] for r in self.records.values())):
            raise ValueError("Invalid original file totals")
        self.records["MANIFEST.json"] = dict(path="MANIFEST.json", bytes=len(original),
            sha256=hashlib.sha256(original).hexdigest())
        self.locations = {}
        self.archive_members = {}
        for row in self.layout["files"]:
            name = safe_name(row["path"])
            if name in self.locations or not row["copies"]:
                raise ValueError("Repeated or unavailable original file")
            self.locations[name] = row["copies"]
            for copy in row["copies"]:
                delivery = copy["delivery"]
                if delivery not in self.roots:
                    raise ValueError("Unknown delivery")
                safe_name(copy["path"])
                if "member" in copy:
                    if delivery != "data":
                        raise ValueError("Data archives must be in the data delivery")
                    member = safe_name(copy["member"])
                    members = self.archive_members.setdefault(copy["path"], set())
                    if member in members:
                        raise ValueError("Repeated archive location")
                    members.add(member)
        if set(self.locations) != set(self.records):
            raise ValueError("The layout does not cover the exact original file set")
        for name, copies in self.locations.items():
            if name.startswith(("code/", "tests/")) and not any(
                    c["delivery"] == "github" and "member" not in c for c in copies):
                raise ValueError("Original source is missing from GitHub: " + name)

    def open_archives(self, stack):
        archives = {}
        for name, expected in self.archive_members.items():
            archive = stack.enter_context(zipfile.ZipFile(local_file(self.roots["data"], name)))
            infos = archive.infolist()
            actual = [safe_name(info.filename) for info in infos]
            if len(actual) != len(set(actual)) or set(actual) != expected:
                raise ValueError("Archive member set differs from layout: " + name)
            for info in infos:
                if (info.is_dir() or stat.S_ISLNK(info.external_attr >> 16)
                        or info.flag_bits & 1):
                    raise ValueError("Unexpected directory, link or encrypted ZIP member")
                if info.filename.endswith((".zst", ".gz", ".zip")) and info.compress_type != zipfile.ZIP_STORED:
                    raise ValueError("Already-compressed data were recompressed")
            archives[name] = archive
        return archives

    def open_copy(self, copy, archives):
        if "member" in copy:
            return archives[copy["path"]].open(copy["member"])
        return local_file(self.roots[copy["delivery"]], copy["path"]).open("rb")

    def verify(self):
        package_checks = {
            name: verify_manifest(root, allow_extra=name == "github")
            for name, root in self.roots.items()
        }
        github_manifest = json.loads((self.roots["github"] / "MANIFEST.json").read_bytes())
        if max(row["bytes"] for row in github_manifest["files"]) >= 100 * 1024**2:
            raise ValueError("A GitHub file reaches the regular Git size limit")
        data_files = list(self.roots["data"].iterdir())
        if (len(data_files) > 500 or any(not p.is_file() for p in data_files)
                or sum(p.stat().st_size for p in data_files) > 20_000_000_000):
            raise ValueError("Data delivery exceeds the Figshare free-account upload limits")
        copies_checked = 0
        original_frontier = None
        with ExitStack() as stack:
            archives = self.open_archives(stack)
            for name, record in self.records.items():
                if name == "data/protocols/pareto_frontier.json":
                    with self.open_copy(self.locations[name][0], archives) as stream:
                        original_frontier = stream.read()
                for copy in self.locations[name]:
                    with self.open_copy(copy, archives) as stream:
                        size, digest = stream_digest(stream)
                    if (size, digest) != (record["bytes"], record["sha256"]):
                        raise ValueError("Original bytes not preserved: " + name)
                    copies_checked += 1
        first = local_file(self.roots["github"], "data/protocols/pareto_frontier.json").read_bytes()
        second = local_file(self.roots["data"], "pareto_frontier.json").read_bytes()
        if first != second:
            raise ValueError("The protocol catalogues differ")
        presentation = None
        current = json.loads(first)
        if "output_representative_presentation" in current:
            from output_factorisation import verify_presentation_update
            if hashlib.sha256(original_frontier).hexdigest() != current[
                    "output_representative_presentation"]["previous_catalogue_sha256"]:
                raise ValueError("The presentation uses a different original catalogue")
            presentation = verify_presentation_update(json.loads(original_frontier), current)
        elif first != original_frontier:
            raise ValueError("Changed catalogue without a presentation certificate")
        return dict(schema="triorthogonal-delivery-verification-v1", status="pass",
            original_files_including_manifest=len(self.records), copies_checked=copies_checked,
            original_bytes=sum(r["bytes"] for r in self.records.values()),
            original_manifest_sha256=self.layout["original_manifest_sha256"],
            every_original_file_preserved=True, all_original_code_in_github=True,
            inner_compressed_streams_unchanged=True, catalogue=catalogue_summary(first),
            packages=package_checks, output_presentation_verification=presentation,
            classification_recomputed=False)

    def assemble(self, destination):
        destination = Path(destination).absolute()
        if destination.exists() or destination.is_symlink():
            raise ValueError("Destination must not exist")
        resolved = destination.resolve()
        if any(resolved.is_relative_to(root) or root.is_relative_to(resolved)
               for root in self.roots.values()):
            raise ValueError("Destination must be outside both deliveries")
        result = self.verify()
        destination.mkdir(parents=True, exist_ok=False)
        with ExitStack() as stack:
            archives = self.open_archives(stack)
            for name in self.records:
                target = destination / name
                target.parent.mkdir(parents=True, exist_ok=True)
                with self.open_copy(self.locations[name][0], archives) as src, target.open("xb") as dst:
                    shutil.copyfileobj(src, dst, 1024 * 1024)
        result["assembled_release"] = verify_manifest(destination)
        return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("verify", "assemble"))
    parser.add_argument("--github", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--data", "--figshare", dest="data", type=Path,
                        required=True, help="Directory containing the downloaded data delivery")
    parser.add_argument("--destination", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.output and any(args.output.resolve().is_relative_to(p.resolve())
                           for p in (args.github, args.data)):
        parser.error("Write audit output outside the two deliveries")
    if args.command == "assemble" and args.destination is None:
        parser.error("assemble requires --destination")
    delivery = Delivery(args.github, args.data)
    result = delivery.assemble(args.destination) if args.command == "assemble" else delivery.verify()
    rendered = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.write_text(rendered, encoding="utf-8")
    print(rendered, end="")


if __name__ == "__main__":
    main()
