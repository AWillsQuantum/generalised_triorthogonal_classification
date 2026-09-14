import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import subprocess
import tempfile
import unittest
import warnings
import zipfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))
from delivery import Delivery, catalogue_summary, safe_name
from verify_release_integrity import create, verify


class DeliveryTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.original = self.root / "original"
        self.github = self.root / "github"
        self.data = self.root / "data"
        for path in (self.original, self.github, self.data):
            path.mkdir()
        self.put(self.original, "code/example.py", b"# Original source\n")
        self.put(self.original, "README.md", b"Original README\n")
        self.put(self.original, "data/space.zst", b"compressed fixture\x00\xff")
        catalogue = dict(maximum_protocol_length=54, minimum_distance=3,
                         equivalence="CNOT+S", intrinsic_outputs_only=True,
                         pareto_partition=["output_id", "exact d_Z"],
                         pareto_objectives=["n", "S"], protocols=[{}] * 74,
                         outputs=[{"representative_gate": "T1"}] * 62)
        self.put(self.original, "data/protocols/pareto_frontier.json", json.dumps(catalogue).encode())
        create(self.original)
        original_raw = (self.original / "MANIFEST.json").read_bytes()
        names = [r["path"] for r in json.loads(original_raw)["files"]] + ["MANIFEST.json"]
        rows = []
        for name in names:
            copies = []
            if name == "data/space.zst":
                with zipfile.ZipFile(self.data / "spaces.zip", "w") as archive:
                    archive.write(self.original / name, name)
                copies.append(dict(delivery="data", path="spaces.zip", member=name))
            else:
                target = {"README.md": "provenance/ORIGINAL_README.md",
                          "MANIFEST.json": "provenance/ORIGINAL_MANIFEST.json"}.get(name, name)
                self.put(self.github, target, (self.original / name).read_bytes())
                copies.append(dict(delivery="github", path=target))
            if name in ("MANIFEST.json", "data/protocols/pareto_frontier.json"):
                target = "ORIGINAL_MANIFEST.json" if name == "MANIFEST.json" else "pareto_frontier.json"
                self.put(self.data, target, (self.original / name).read_bytes())
                copies.append(dict(delivery="data", path=target))
            rows.append(dict(path=name, copies=copies))
        self.layout = dict(schema="triorthogonal-two-delivery-layout-v2", files=rows,
                           original_manifest_sha256=hashlib.sha256(original_raw).hexdigest())
        self.write_layout()

    @staticmethod
    def put(root, name, raw):
        target = root / name
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(raw)

    def manifests(self):
        create(self.github)
        create(self.data)

    def write_layout(self):
        for root in (self.github, self.data):
            (root / "delivery_layout.json").write_text(json.dumps(self.layout), encoding="ascii")
        self.manifests()

    def delivery(self):
        return Delivery(self.github, self.data)

    def test_preservation_and_exact_reassembly(self):
        result = self.delivery().verify()
        self.assertTrue(result["every_original_file_preserved"])
        destination = self.root / "assembled"
        self.delivery().assemble(destination)
        self.assertEqual(verify(destination), verify(self.original))
        self.assertEqual({p.relative_to(destination) for p in destination.rglob("*") if p.is_file()},
                         {p.relative_to(self.original) for p in self.original.rglob("*") if p.is_file()})

    def test_missing_original_mapping(self):
        self.layout["files"].pop()
        self.write_layout()
        with self.assertRaisesRegex(ValueError, "exact original file set"):
            self.delivery()

    def test_changed_catalogue_despite_new_package_hashes(self):
        self.put(self.data, "pareto_frontier.json", b"{}")
        self.manifests()
        with self.assertRaisesRegex(ValueError, "not preserved"):
            self.delivery().verify()

    def test_changed_archive_despite_new_package_hashes(self):
        with zipfile.ZipFile(self.data / "spaces.zip", "w") as archive:
            archive.writestr("data/space.zst", b"different")
        self.manifests()
        with self.assertRaisesRegex(ValueError, "not preserved"):
            self.delivery().verify()

    def test_unsafe_paths(self):
        for name in ("", ".", "../outside", "/absolute", "C:/absolute", "a\\b", "./local", "a//b", "a\x00b"):
            with self.subTest(name=name), self.assertRaises(ValueError):
                safe_name(name)

    def test_unexpected_archive_member(self):
        with zipfile.ZipFile(self.data / "spaces.zip", "a") as archive:
            archive.writestr("../outside", b"unsafe")
        self.manifests()
        with self.assertRaisesRegex(ValueError, "Unsafe relative path"):
            self.delivery().verify()

    def test_duplicate_archive_member(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            with zipfile.ZipFile(self.data / "spaces.zip", "a") as archive:
                archive.writestr("data/space.zst", b"duplicate")
        self.manifests()
        with self.assertRaisesRegex(ValueError, "member set"):
            self.delivery().verify()

    def test_no_recompression_of_inner_data(self):
        with zipfile.ZipFile(self.data / "spaces.zip", "w", compression=zipfile.ZIP_DEFLATED) as archive:
            archive.write(self.original / "data/space.zst", "data/space.zst")
        self.manifests()
        with self.assertRaisesRegex(ValueError, "recompressed"):
            self.delivery().verify()

    def test_destination_inside_delivery_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "outside both"):
            self.delivery().assemble(self.github / "output")

    def test_existing_destination_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "must not exist"):
            self.delivery().assemble(self.original)

    def test_bad_scope_is_rejected(self):
        raw = json.loads((self.data / "pareto_frontier.json").read_bytes())
        raw["equivalence"] = "full Clifford"
        with self.assertRaisesRegex(ValueError, "scope or counts"):
            catalogue_summary(json.dumps(raw).encode())

    def test_mismatched_layouts(self):
        (self.data / "delivery_layout.json").write_text("{}")
        with self.assertRaisesRegex(ValueError, "layouts differ"):
            self.delivery()

    def test_data_cli_aliases(self):
        for option in ("--data", "--figshare"):
            with self.subTest(option=option):
                result = subprocess.run(
                    [sys.executable, "-B", str(ROOT / "code/delivery.py"), "verify",
                     "--github", str(self.github), option, str(self.data)],
                    capture_output=True, text=True)
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertTrue(json.loads(result.stdout)["every_original_file_preserved"])

    def test_builder_excludes_author_publication_instructions(self):
        self.put(self.original, "code/verify_release_integrity.py",
                 (ROOT / "code/verify_release_integrity.py").read_bytes())
        create(self.original)
        spec = importlib.util.spec_from_file_location(
            "delivery_builder", ROOT / "packaging/build_deliveries.py")
        builder = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(builder)
        destination = self.root / "rebuilt"
        builder.build(self.original, destination, ROOT)
        self.assertEqual((destination / "github/CITATION.cff").read_bytes(),
                         (ROOT / "CITATION.cff").read_bytes())
        self.assertEqual((destination / "figshare/README.md").read_bytes(),
                         (ROOT / "packaging/figshare_README.md").read_bytes())
        self.assertIn("Adam Wills is the citation author",
                      (destination / "figshare/README.md").read_text(encoding="utf-8"))
        self.assertIn("ChatGPT (OpenAI)",
                      (destination / "figshare/README.md").read_text(encoding="utf-8"))
        for delivery_name in ("github", "figshare"):
            self.assertEqual((destination / delivery_name / "LICENSE_DATA.md").read_bytes(),
                             (ROOT / "LICENSE_DATA.md").read_bytes())
        for path in destination.rglob("*"):
            self.assertNotIn(path.name, ("PUBLISHING.md", "LICENSING.md",
                                         "FIGSHARE_UPLOAD.md"))
            if path.name.endswith("README.md"):
                text = path.read_text(encoding="utf-8")
                self.assertNotIn("PUBLISHING.md", text)
                self.assertNotIn("Upload all top-level files", text)
        self.assertTrue(Delivery(destination / "github", destination / "figshare")
                        .verify()["every_original_file_preserved"])


if __name__ == "__main__":
    unittest.main()
