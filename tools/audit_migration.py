#!/usr/bin/env python3
"""Verify copied source, recorded hashes, and unchanged original source files."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess

from prepare_nwchem import LINK_BLOCK, LINK_MARKER

parser = argparse.ArgumentParser(description=__doc__)
project = Path(__file__).resolve().parents[1]
parser.add_argument("--manifest", type=Path, default=project / "migration_manifest.json")
parser.add_argument("--destination", type=Path, default=project / "upstream")
args = parser.parse_args()
manifest = json.loads(args.manifest.read_text())
source = Path(manifest["source"])
changes = {item["file"]: item for item in manifest["changed_files"]}
raw = subprocess.check_output(["git", "-C", str(source), "ls-files", "-z", "src"])
checked = 0
for encoded in raw.split(b"\0"):
    if not encoded:
        continue
    relative = os.fsdecode(encoded)
    original, copied = source / relative, args.destination / relative
    if original.is_symlink():
        assert copied.is_symlink() and os.readlink(original) == os.readlink(copied), relative
        continue
    assert original.is_file() and copied.is_file(), relative
    original_data, copied_data = original.read_bytes(), copied.read_bytes()
    if relative in changes:
        item = changes[relative]
        assert hashlib.sha256(original_data).hexdigest() == item["original_sha256"], "Original changed: " + relative
        assert hashlib.sha256(copied_data).hexdigest() == item["ported_sha256"], "Port changed: " + relative
        assert [len(x) for x in original_data.splitlines()] == [len(x) for x in copied_data.splitlines()], "Fixed-form columns changed: " + relative
    elif relative == "src/config/makefile.h":
        expected = original_data.decode().replace(LINK_MARKER, LINK_BLOCK.rstrip()).encode()
        assert copied_data == expected, "Integration link configuration differs"
    else:
        assert original_data == copied_data, "Unexpected source change: " + relative
    checked += 1
result = {"source_files_audited": checked, "changed_fortran_files_verified": len(changes),
          "original_files_preserved": True, "fixed_form_line_lengths_preserved": True}
print(json.dumps(result, indent=2))
