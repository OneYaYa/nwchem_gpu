import importlib.util
import contextlib
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

spec = importlib.util.spec_from_file_location("prepare", Path(__file__).parents[1] / "tools/prepare_nwchem.py")
prepare = importlib.util.module_from_spec(spec)
spec.loader.exec_module(prepare)


class RewriteTest(unittest.TestCase):
    def test_fixed_form_and_comments(self):
        original = ("c     call dgemm('n','n')\n"
                    "      CALL DGEMM('n','n',m,n,k,1d0,a,m,b,k,0d0,c,m) ! call dgemm\n"
                    "     & call dgemm('t','n')\n"
                    "      call my_dgemm(a)\n"
                    "#define NO_CHANGE call dgemm\n")
        result, count = prepare.rewrite_fortran(original)
        self.assertEqual(count, 2)
        self.assertEqual([len(x) for x in result.splitlines()], [len(x) for x in original.splitlines()])
        self.assertIn("c     call dgemm", result)
        self.assertIn("CALL nwgmm", result)
        self.assertIn("! call dgemm", result)
        self.assertIn("call my_dgemm", result)
        self.assertEqual(prepare.rewrite_fortran(result)[1], 0)

    def test_recursive_include(self):
        text = ('      Subroutine xc_dirac(a)\n      Subroutine xc_dirac_d2(a)\n'
                '      Subroutine xc_dirac_d3(a)\n#include "xc_dirac.F"\n')
        result = prepare.port_dirac(text)
        self.assertIn('include "xc_dirac.F"', result)
        self.assertIn('xc_origx_d3', result)
        self.assertEqual(len(text), len(result))

    def test_snapshot_and_source_preservation(self):
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            source, destination = root / "original", root / "work" / "upstream"
            files = {
                "src/config/makefile.h": prepare.LINK_MARKER + "\n",
                "src/nwdft/xc/xc_dirac.F": '      subroutine xc_dirac(a)\n'
                '      subroutine xc_dirac_d2(a)\n      subroutine xc_dirac_d3(a)\n#include "xc_dirac.F"\n',
                "src/tce/a.F": "      call DGEMM('n','t')\n",
                "src/md/preserve.F": "      call unrelated(a)\n",
                "QA/reference.out": "saved scientific output\n",
            }
            for name, value in files.items():
                target = source / name
                target.parent.mkdir(parents=True, exist_ok=True)
                target.write_text(value)
            library = root / "libnwchem_gpu.so"
            library.write_text("test fixture")

            def git_output(command, **kwargs):
                if "ls-files" in command:
                    return ("\0".join(files)+"\0").encode()
                return "test-revision\n"

            with patch.object(prepare.subprocess, "check_output", side_effect=git_output):
                with contextlib.redirect_stdout(io.StringIO()):
                    prepare.prepare(source, destination, library)
            for name, value in files.items():
                self.assertEqual((source / name).read_text(), value)
            self.assertIn("call nwgmm", (destination / "src/tce/a.F").read_text())
            self.assertEqual((destination / "src/md/preserve.F").read_text(), files["src/md/preserve.F"])
            config = (destination / "src/config/makefile.h").read_text()
            self.assertIn("$(CORE_LIBS) $(NWGPU_LIBRARY)", config)
            self.assertIn("-Wl,-rpath", config)
            manifest = json.loads((destination.parent / "migration_manifest.json").read_text())
            self.assertEqual(manifest["modules"]["tce"]["gemm_calls_routed"], 1)
            self.assertEqual(manifest["modules"]["md"]["gemm_calls_routed"], 0)
            self.assertTrue((destination / "QA/reference.out").exists())
            with self.assertRaisesRegex(ValueError, "already exists"):
                prepare.prepare(source, destination, library)


if __name__ == "__main__":
    unittest.main()
