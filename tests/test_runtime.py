from __future__ import annotations

import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

APP = Path(__file__).resolve().parents[1] / "app"
sys.path.insert(0, str(APP))

import geometry_native
import scatmech_paths as paths


class RuntimeTests(unittest.TestCase):
    def setUp(self):
        self.context = tempfile.TemporaryDirectory(prefix="scatmech-runtime-")
        self.root = Path(self.context.name)
        self.addCleanup(self.context.cleanup)
        self.environment = patch.dict(os.environ, {
            "SCATMECH_BIN": "", "SCATMECH_DATA_DIR": "", "SCATMECH_CACHE_DIR": "",
            "SCATMECH": "", "LOCALAPPDATA": str(self.root / "User Data"),
        })
        self.environment.start()
        self.addCleanup(self.environment.stop)

    def test_frozen_resource_and_writable_data_locations_are_separate(self):
        bundle = self.root / "read-only-install" / "_internal"
        with patch.object(sys, "frozen", True, create=True), \
                patch.object(sys, "_MEIPASS", str(bundle), create=True), \
                patch.object(sys, "platform", "win32"):
            self.assertEqual(paths.get_scatmech_bin(), bundle / "solvers")
            self.assertEqual(paths.get_data_dir(), self.root / "User Data" / "SCATMECH-GUI" / "DATA")
            self.assertEqual(paths.get_cache_dir(), self.root / "User Data" / "SCATMECH-GUI" / "cache")

    def test_explicit_paths_override_packaged_defaults(self):
        with patch.object(sys, "frozen", True, create=True), patch.dict(os.environ, {
            "SCATMECH_BIN": str(self.root / "native"),
            "SCATMECH_DATA_DIR": str(self.root / "results"),
            "SCATMECH_CACHE_DIR": str(self.root / "cache"),
        }):
            self.assertEqual(paths.get_scatmech_bin(), (self.root / "native").resolve())
            self.assertEqual(paths.get_data_dir(), (self.root / "results").resolve())
            self.assertEqual(paths.get_cache_dir(), (self.root / "cache").resolve())

    def test_frozen_launch_uses_writable_workdir_and_solver_search_path(self):
        previous = Path.cwd()
        self.addCleanup(os.chdir, previous)
        bundle = self.root / "install" / "_internal"
        with patch.object(sys, "frozen", True, create=True), \
                patch.object(sys, "_MEIPASS", str(bundle), create=True), \
                patch.object(sys, "platform", "win32"):
            paths.configure_runtime(change_directory=True)
            self.assertEqual(Path.cwd(), paths.get_data_dir().resolve())
            self.assertIn(str(bundle / "solvers"), os.environ["SCATMECH"].split(os.pathsep))
            self.assertFalse(bundle.exists())

    def test_frozen_worker_uses_sibling_console_executable(self):
        name = "SCATMECH-Worker.exe" if os.name == "nt" else "SCATMECH-Worker"
        worker = self.root / name
        worker.touch()
        with patch.object(sys, "frozen", True, create=True), \
                patch.object(sys, "executable", str(self.root / "SCATMECH-GUI.exe")):
            program, args = paths.worker_command(APP / "custom_profile_worker.py", "request.json", "output folder")
            self.assertEqual(program, str(worker.resolve()))
            self.assertEqual(args, ["--custom-profile-worker", "request.json", "output folder"])
            worker.unlink()
            with self.assertRaisesRegex(FileNotFoundError, "complete application folder"):
                paths.worker_command(APP / "geometry_worker.py", "request", "output")

    def test_source_worker_uses_python_and_script(self):
        with patch.object(sys, "frozen", False, create=True):
            program, args = paths.worker_command(APP / "geometry_worker.py", "in.json", "out.npz")
        expected = Path(sys.executable)
        if expected.name.lower() == "pythonw.exe" and expected.with_name("python.exe").is_file():
            expected = expected.with_name("python.exe")
        self.assertEqual(program, str(expected))
        self.assertEqual(args, [str(APP / "geometry_worker.py"), "in.json", "out.npz"])

    def test_prebuilt_geometry_exporter_needs_no_compiler_or_sources(self):
        helper = self.root / ("geometry_export.exe" if os.name == "nt" else "geometry_export")
        helper.touch()
        with patch.object(geometry_native, "get_scatmech_bin", return_value=self.root), \
                patch.object(geometry_native, "_source_dir") as sources, \
                patch.object(geometry_native, "_compiler") as compiler:
            self.assertEqual(geometry_native.ensure_exporter(), helper)
            sources.assert_not_called()
            compiler.assert_not_called()

    def test_missing_frozen_helper_does_not_attempt_a_build(self):
        with patch.object(geometry_native, "get_scatmech_bin", return_value=self.root), \
                patch.object(geometry_native, "is_frozen", return_value=True), \
                patch.object(geometry_native, "_source_dir") as sources:
            with self.assertRaisesRegex(RuntimeError, "packaged geometry helper is missing"):
                geometry_native.ensure_exporter()
            sources.assert_not_called()

    def test_custom_profile_cache_uses_windows_alias_for_spaced_user_path(self):
        with patch.object(sys, "platform", "win32"), \
                patch.object(paths, "native_path", return_value=self.root / "short-cache") as shorten:
            self.assertEqual(paths.get_native_cache_dir(), self.root / "short-cache")
            directory = shorten.call_args.args[0]
            self.assertTrue(directory.is_dir())
            self.assertIn("User Data", str(directory))

    def test_native_cache_reports_override_when_short_paths_are_unavailable(self):
        with patch.object(sys, "platform", "win32"), \
                patch.object(paths, "native_path", side_effect=lambda path: path):
            with self.assertRaisesRegex(RuntimeError, "SCATMECH_CACHE_DIR"):
                paths.get_native_cache_dir()

    def test_worker_entry_rejects_invalid_command_without_starting_gui(self):
        result = subprocess.run([sys.executable, str(APP / "worker_entry.py"), "--invalid"],
                                capture_output=True, text=True, check=False)
        self.assertEqual(result.returncode, 2)
        self.assertIn("Usage: SCATMECH-Worker", result.stderr)


if __name__ == "__main__":
    unittest.main()
