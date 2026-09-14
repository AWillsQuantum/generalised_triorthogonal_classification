"""Create or verify the complete file manifest; this is an integrity check, not a proof."""

import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath


MANIFEST = "MANIFEST.json"
LIMIT = 50_000_000_000


def files(root):
    result = {}
    for path in sorted(root.rglob("*")):
        if path.is_symlink():
            raise ValueError("Symbolic link in release: "+str(path))
        if path.is_file():
            name = path.relative_to(root).as_posix()
            if name != MANIFEST:
                result[name] = path
    return result


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def create(root):
    root = Path(root).resolve()
    inventory = files(root)
    if any("__pycache__" in PurePosixPath(name).parts or name.endswith(".pyc") for name in inventory):
        raise ValueError("Remove generated Python bytecode before publishing the release")
    records = [dict(path=name, bytes=path.stat().st_size, sha256=digest(path))
               for name, path in inventory.items()]
    result = dict(schema="triorthogonal-release-file-manifest-v1", algorithm="sha256",
                  excluded_files=[MANIFEST], file_count=len(records),
                  payload_bytes=sum(row["bytes"] for row in records), files=records)
    content = (json.dumps(result, indent=2)+"\n").encode("ascii")
    if result["payload_bytes"]+len(content) >= LIMIT:
        raise ValueError("The distribution must be smaller than 50 GB")
    (root / MANIFEST).write_bytes(content)
    return result


def verify(root, allow_extra=False):
    root = Path(root).resolve()
    raw = (root / MANIFEST).read_bytes()
    data = json.loads(raw)
    if (data.get("schema") != "triorthogonal-release-file-manifest-v1"
            or data.get("algorithm") != "sha256" or data.get("excluded_files") != [MANIFEST]):
        raise ValueError("Unknown release manifest")
    records = data["files"]
    declared = {}
    for row in records:
        name = row["path"]
        if not isinstance(name, str):
            raise ValueError("Invalid manifest path")
        path = PurePosixPath(name)
        if (not name or path.is_absolute() or str(path) != name or ".." in path.parts
                or ":" in name or "\\" in name or name == MANIFEST or name in declared):
            raise ValueError("Unsafe or repeated manifest path")
        if type(row["bytes"]) is not int or row["bytes"] < 0:
            raise ValueError("Invalid file size")
        if (not isinstance(row["sha256"], str) or len(row["sha256"]) != 64
                or any(ch not in "0123456789abcdef" for ch in row["sha256"])):
            raise ValueError("Invalid SHA-256 digest")
        declared[name] = row
    actual = files(root)
    missing, extra = set(declared)-set(actual), set(actual)-set(declared)
    if missing or extra and not allow_extra:
        raise ValueError(f"File set mismatch: missing={sorted(missing)}, extra={sorted(extra)}")
    for name, row in declared.items():
        if actual[name].stat().st_size != row["bytes"] or digest(actual[name]) != row["sha256"]:
            raise ValueError("Changed release file: "+name)
    size = sum(row["bytes"] for row in records)
    if data["file_count"] != len(records) or data["payload_bytes"] != size or size+len(raw) >= LIMIT:
        raise ValueError("Invalid file count or release size")
    return dict(schema="triorthogonal-release-integrity-verification-v1", status="pass",
                files=len(records), payload_bytes=size, distribution_bytes=size+len(raw),
                manifest_sha256=hashlib.sha256(raw).hexdigest(), extra_files=len(extra),
                classification_completeness_established=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--create", action="store_true", help="Replace the manifest with hashes of current files")
    parser.add_argument("--allow-extra", action="store_true", help="Permit additional local build or audit files")
    args = parser.parse_args()
    if args.create:
        create(args.root)
    print(json.dumps(verify(args.root, args.allow_extra), indent=2))


if __name__ == "__main__":
    main()
