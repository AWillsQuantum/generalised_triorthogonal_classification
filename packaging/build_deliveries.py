"""Build code and data deliveries from a verified combined classification release."""

import argparse
from contextlib import ExitStack
import hashlib
import json
from pathlib import Path
import shutil
import sys
import zipfile


EXTRAS = ("README.md", "CITATION.cff", "LICENSE_DATA.md", ".gitignore", ".gitattributes",
          "code/delivery.py", "delivery_tests/test_delivery.py",
          "packaging/figshare_README.md", "packaging/figshare_DATA_FORMAT.md")


def copy_file(source, target):
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(source, target)


def build(source, destination, templates):
    source, destination, templates = (Path(p).resolve() for p in (source, destination, templates))
    if destination.exists() or destination.is_relative_to(source) or source.is_relative_to(destination):
        raise ValueError("Use a new staging directory outside the original release")
    sys.path.insert(0, str(source / "code"))
    from verify_release_integrity import create, verify
    baseline = verify(source)
    raw = (source / "MANIFEST.json").read_bytes()
    original = json.loads(raw)
    records = original["files"] + [dict(path="MANIFEST.json", bytes=len(raw),
                                       sha256=hashlib.sha256(raw).hexdigest())]
    github, data = destination / "github", destination / "figshare"
    github.mkdir(parents=True)
    data.mkdir()
    rows = []
    with ExitStack() as stack:
        archives = {name: stack.enter_context(zipfile.ZipFile(data / name, "w", allowZip64=True))
                    for name in ("space_catalogues.zip", "classification_evidence.zip")}
        for record in records:
            name = record["path"]
            copies = []
            if name.startswith("data/spaces/"):
                archive = "space_catalogues.zip"
            elif name.startswith(("certificates/", "data/")) and not name.startswith(("data/protocols/", "data/examples/")):
                archive = "classification_evidence.zip"
            else:
                archive = None
            if archive:
                info = zipfile.ZipInfo(name, date_time=(2026, 9, 12, 0, 0, 0))
                info.external_attr = 0o100644 << 16
                info.compress_type = zipfile.ZIP_STORED if name.endswith((".zst", ".zip", ".gz")) else zipfile.ZIP_DEFLATED
                with (source / name).open("rb") as src, archives[archive].open(info, "w", force_zip64=True) as dst:
                    shutil.copyfileobj(src, dst, 1024 * 1024)
                copies.append(dict(delivery="data", path=archive, member=name))
            if not archive or name in ("data/spaces/index.json", "data/spaces/c16/m04/spaces_000000000.utspace.zst"):
                target = {"README.md": "provenance/ORIGINAL_README.md",
                          "MANIFEST.json": "provenance/ORIGINAL_MANIFEST.json"}.get(name, name)
                copy_file(source / name, github / target)
                copies.append(dict(delivery="github", path=target))
            if name in ("data/protocols/pareto_frontier.json", "MANIFEST.json"):
                target = "pareto_frontier.json" if name.startswith("data/") else "ORIGINAL_MANIFEST.json"
                copy_file(source / name, data / target)
                copies.append(dict(delivery="data", path=target))
            rows.append(dict(path=name, copies=copies))
    layout = dict(schema="triorthogonal-two-delivery-layout-v2",
                  original_manifest_sha256=baseline["manifest_sha256"], files=rows)
    layout_raw = (json.dumps(layout, indent=2) + "\n").encode("ascii")
    for root in (github, data):
        (root / "delivery_layout.json").write_bytes(layout_raw)
    for name in EXTRAS:
        copy_file(templates / name, github / name)
    copy_file(Path(__file__), github / "packaging/build_deliveries.py")
    for name in ("README", "DATA_FORMAT"):
        copy_file(templates / f"packaging/figshare_{name}.md", data / f"{name}.md")
    copy_file(templates / "LICENSE_DATA.md", data / "LICENSE_DATA.md")
    for root in (github, data):
        create(root)
    return baseline


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--destination", type=Path, required=True)
    parser.add_argument("--templates", type=Path, required=True,
                        help="GitHub delivery containing the documentation and packaging helpers")
    args = parser.parse_args()
    print(json.dumps(build(args.source, args.destination, args.templates), indent=2))


if __name__ == "__main__":
    main()
