#!/usr/bin/env python3
"""Create a separate, complete tracked NWChem tree with explicit GPU call sites.

Never edits the source tree, overwrites a destination, or changes BLAS symbols.
The source snapshot contains every tracked module, QA fixture, and build file.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import shlex
import subprocess

FORTRAN = {".f", ".f90", ".f95", ".fh", ".finc"}
CALL = re.compile(r"\bcall[ \t]+dgemm\b", re.IGNORECASE)
DIRAC = re.compile(r"(\bsubroutine[ \t]+)xc_dirac(_d[23])?\b", re.IGNORECASE)
LINK_MARKER = "LIBS = $(NW_MODULE_LIBS) $(CORE_LIBS)"
LINK_BLOCK = """# C++/CUDA backend; resolve calls after module and core archives.
ifndef NWGPU_LIBRARY
  $(error Set NWGPU_LIBRARY to the absolute path of libnwchem_gpu.so or .dylib)
endif
LDFLAGS += -Wl,-rpath,$(dir $(NWGPU_LIBRARY))
LIBS = $(NW_MODULE_LIBS) $(CORE_LIBS) $(NWGPU_LIBRARY)
"""


def rewrite_fortran(text):
    count = 0
    output = []
    for line in text.splitlines(keepends=True):
        # Preserve fixed-form comments, inline comments, and every line length.
        if line and line[0] in "cC*!#":
            output.append(line)
            continue
        code, sep, comment = line.partition("!")
        code, changed = CALL.subn(lambda m: re.sub("dgemm$", "nwgmm", m[0], flags=re.I), code)
        count += changed
        output.append(code + sep + comment)
    return "".join(output), count


def port_dirac(text):
    # Keep the original derivative implementations compiled under private names;
    # the recursive #include must still refer to the original filename.
    renamed, count = DIRAC.subn(lambda m: m[1] + "xc_origx" + (m[2] or ""), text)
    if count != 3:
        raise ValueError("Expected exactly three xc_dirac routine declarations")
    return renamed


def prepare(source, destination, library):
    source = source.resolve()
    destination = destination.resolve()
    library = library.resolve()
    if destination.exists():
        raise ValueError("Destination already exists; choose a new empty path")
    if not library.is_file():
        raise ValueError("Build the backend first; library does not exist: " + str(library))
    # GNU make cannot safely represent spaces/metacharacters in raw link paths.
    if not re.fullmatch(r"[A-Za-z0-9_./+:-]+", str(library)):
        raise ValueError("Backend library path must contain only letters, digits, _, ., /, +, :, -")
    if destination == source or destination in source.parents:
        raise ValueError("Destination cannot contain the source")
    config = source / "src/config/makefile.h"
    if config.read_text().count(LINK_MARKER) != 1:
        raise ValueError("Unsupported upstream makefile layout")
    # Use Git's index, including working-tree edits of tracked files; exclude build artifacts.
    raw = subprocess.check_output(["git", "-C", str(source), "ls-files", "-z"])
    paths = [Path(os.fsdecode(p)) for p in raw.split(b"\0") if p]
    destination.mkdir(parents=True)
    modules = {}
    changes = []
    for relative in paths:
        original = source / relative
        if not original.is_file() and not original.is_symlink():
            continue
        # A later commit may track this project; never recursively include it.
        if relative.parts[0] == "nwchem_gpu":
            continue
        target = destination / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        if original.is_symlink():
            target.symlink_to(os.readlink(original))
            continue
        shutil.copy2(original, target)
        if relative.parts[0] != "src":
            continue
        module = relative.parts[1] if len(relative.parts) > 2 else "(build)"
        info = modules.setdefault(module, {"files": 0, "fortran_files": 0, "gemm_calls_routed": 0})
        info["files"] += 1
        if relative.suffix.lower() not in FORTRAN:
            continue
        info["fortran_files"] += 1
        data = original.read_bytes()
        text = data.decode("latin-1")
        result, count = rewrite_fortran(text)
        if str(relative) == "src/nwdft/xc/xc_dirac.F":
            result = port_dirac(result)
        if result != text:
            target.write_bytes(result.encode("latin-1"))
            info["gemm_calls_routed"] += count
            changes.append({"file": str(relative), "gemm_calls": count,
                            "original_sha256": hashlib.sha256(data).hexdigest(),
                            "ported_sha256": hashlib.sha256(target.read_bytes()).hexdigest()})
    config_target = destination / "src/config/makefile.h"
    config_target.write_text(config_target.read_text().replace(LINK_MARKER, LINK_BLOCK.rstrip()))
    revision = subprocess.check_output(["git", "-C", str(source), "rev-parse", "HEAD"], text=True).strip()
    config_path = library.parent / "nwchem_gpu_config.json"
    backend_config = json.loads(config_path.read_text()) if config_path.is_file() else None
    manifest = {"source": str(source), "source_revision": revision,
                "source_snapshot": "Tracked working-tree files, including local tracked edits; no build artifacts or untracked tutorials",
                "backend_library": str(library), "backend_configuration": backend_config,
                "modules": dict(sorted(modules.items())),
                "changed_files": changes,
                "limitations": ["Fortran remains for all routines except xc_dirac/_d2/_d3 and the routed DGEMM calls",
                                "Indirect BLAS calls, external libraries and runtime-generated source are not rewritten",
                                "Full NWChem build, QA and GPU speedup must be validated separately"]}
    (destination.parent / "migration_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    (destination / "nwchem_gpu_env.sh").write_text(
        "# Source this after selecting the normal NWChem build environment.\n"
        + "export NWCHEM_TOP=" + shlex.quote(str(destination)) + "\n"
        + "export NWGPU_LIBRARY=" + shlex.quote(str(library)) + "\n")
    print(json.dumps({"destination": str(destination), "source_groups": len(modules),
                      "changed_files": len(changes),
                      "gemm_calls_routed": sum(m["gemm_calls_routed"] for m in modules.values())}, indent=2))


def main():
    project = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=project.parent)
    parser.add_argument("--destination", type=Path, default=project / "upstream")
    parser.add_argument("--library", required=True, type=Path)
    args = parser.parse_args()
    prepare(args.source, args.destination, args.library)


if __name__ == "__main__":
    main()
