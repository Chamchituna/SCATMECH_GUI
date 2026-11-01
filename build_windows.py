from __future__ import annotations

import importlib.metadata
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import zipfile

ROOT = Path(__file__).resolve().parent
BUILD = ROOT / "build"
DIST = ROOT / "dist"
NAMES = ("brdfprog", "mieprog", "reflectprog", "rcwprog", "geometry_export")


def run(command, *, timeout=None, env=None):
    print("\n> " + subprocess.list2cmdline(list(map(str, command))), flush=True)
    with (BUILD / "build.log").open("a", encoding="utf-8") as log:
        log.write("\n> " + subprocess.list2cmdline(list(map(str, command))) + "\n")
        log.flush()

        if timeout is not None:
            result = subprocess.run(command, cwd=ROOT, env=env, capture_output=True,
                                    text=True, errors="replace", timeout=timeout)
            print(result.stdout, end="", flush=True)
            print(result.stderr, end="", flush=True)
            log.write(result.stdout + result.stderr)
            if result.returncode:
                raise RuntimeError(f"Command failed with exit code {result.returncode}")
        else:
            with subprocess.Popen(command, cwd=ROOT, env=env, stdout=subprocess.PIPE,
                                  stderr=subprocess.STDOUT, text=True, errors="replace") as process:
                for line in process.stdout:
                    print(line, end="", flush=True)
                    log.write(line)
                if process.wait():
                    raise RuntimeError("Command failed; see build/build.log")


def collect_notices(destination):
    destination.mkdir(parents=True, exist_ok=True)
    for filename in ("LICENSE", "THIRD_PARTY_NOTICES.md"):
        shutil.copy2(ROOT / filename, destination / filename)
    for filename in ("SCATMECH_LICENSE.txt", "THIRD_PARTY_NOTICES.txt", "SOURCE_NOTICES.txt"):
        shutil.copy2(ROOT / "native_src" / filename, destination / filename)
    inventory = []
    for dist in sorted(importlib.metadata.distributions(), key=lambda d: d.metadata["Name"].lower()):
        name = dist.metadata["Name"]
        inventory.append({"name": name, "version": dist.version,
                          "license": dist.metadata.get("License-Expression") or dist.metadata.get("License"),
                          "homepage": dist.metadata.get("Home-page")})
        for entry in dist.files or []:
            if any(word in entry.name.lower() for word in ("license", "copying", "copyright", "notice")):
                source = Path(dist.locate_file(entry))
                if source.is_file() and source.suffix not in {".py", ".pyc", ".pyd", ".dll"}:
                    target = destination / name / str(entry).replace("/", "_").replace("\\", "_")
                    target.parent.mkdir(parents=True, exist_ok=True)
                    shutil.copy2(source, target)
    (destination / "dependency-inventory.json").write_text(json.dumps(inventory, indent=2), encoding="utf-8")
    python_license = Path(sys.base_prefix) / "LICENSE.txt"
    if python_license.is_file():
        shutil.copy2(python_license, destination / "PYTHON_LICENSE.txt")


def source_archive(target):

    allowed_roots = {"app", "native_src", "tests"}
    allowed_files = {"BUILD_WINDOWS.cmd", "build_windows.py", "SCATMECH-GUI.spec",
                     "requirements-build.txt", "README.md", "LICENSE", "THIRD_PARTY_NOTICES.md",
                     "SOURCE_MANIFEST.json", "VALIDATION.md", ".gitignore",
                     "BUILD_WINDOWS.md", "requirements.txt"}
    with zipfile.ZipFile(target, "w", zipfile.ZIP_DEFLATED) as archive:
        for path in sorted(ROOT.rglob("*")):
            rel = path.relative_to(ROOT)
            if not path.is_file() or not (rel.parts[0] in allowed_roots or str(rel) in allowed_files):
                continue
            if any(part in {"__pycache__", "solvers", "DATA", ".runner"} for part in rel.parts):
                continue
            if path.suffix in {".pyc", ".o", ".obj", ".exe", ".dll"} or path.name == ".DS_Store":
                continue
            archive.write(path, Path("windows_distribution") / rel)


def main():
    if sys.platform != "win32":
        raise RuntimeError("Run BUILD_WINDOWS.cmd on Windows 10/11 x64. Windows executables must be built on Windows.")
    if sys.version_info[:2] != (3, 12) or platform.architecture()[0] != "64bit" or platform.machine().lower() not in {"amd64", "x86_64"}:
        raise RuntimeError("Install CPython 3.12 for Windows x64 (not 32-bit or ARM64), then run BUILD_WINDOWS.cmd.")
    os.chdir(ROOT)
    BUILD.mkdir(exist_ok=True)
    venv_python = ROOT / ".build-venv" / "Scripts" / "python.exe"
    if "--inside-venv" not in sys.argv:
        (BUILD / "build.log").write_text("SCATMECH GUI Windows build\n", encoding="utf-8")
        if not venv_python.is_file():
            run([sys.executable, "-m", "venv", str(ROOT / ".build-venv")])
        run([str(venv_python), "-m", "pip", "install", "--upgrade", "pip"])
        run([str(venv_python), "-m", "pip", "install", "--only-binary=:all:", "-r", "requirements-build.txt"])
        return subprocess.call([str(venv_python), str(Path(__file__)), "--inside-venv"], cwd=ROOT)

    cmake = Path(sys.executable).parent / "cmake.exe"
    if not cmake.is_file():
        raise RuntimeError("CMake is missing from .build-venv; rerun BUILD_WINDOWS.cmd.")
    run([str(cmake), "-S", "native_src", "-B", "build/native", "-G", "Visual Studio 17 2022", "-A", "x64"])
    run([str(cmake), "--build", "build/native", "--config", "Release", "--parallel", str(min(8, os.cpu_count() or 2))])
    solvers = ROOT / "app" / "solvers"
    solvers.mkdir(exist_ok=True)
    for name in NAMES:
        shutil.copy2(BUILD / "native" / "bin" / (name + ".exe"), solvers / (name + ".exe"))
    run([sys.executable, "-m", "PyInstaller", "--noconfirm", "--clean", "SCATMECH-GUI.spec"])
    package = DIST / "SCATMECH-GUI"
    report = BUILD / "smoke-test.json"
    report.unlink(missing_ok=True)
    environment = os.environ.copy()
    for name in ("SCATMECH_BIN", "SCATMECH_SOURCE", "SCATMECH_SOURCE_DIR", "SCATMECH", "SCATMECH_CACHE_DIR"):
        environment.pop(name, None)
    environment["SCATMECH_DATA_DIR"] = str(BUILD / "smoke-data")
    run([str(package / "SCATMECH-GUI.exe"), "--smoke-test", str(report)], timeout=240, env=environment)
    if not report.is_file() or not json.loads(report.read_text(encoding="utf-8")).get("passed"):
        raise RuntimeError("Packaged smoke test failed; inspect build/smoke-test.json and build/build.log.")
    shutil.copy2(ROOT / "README.md", package / "README.md")
    shutil.copy2(ROOT / "BUILD_WINDOWS.md", package / "BUILD_WINDOWS.md")
    shutil.copy2(report, package / "BUILD_VERIFICATION.json")
    collect_notices(package / "licenses")
    resolved = subprocess.check_output([sys.executable, "-m", "pip", "freeze"], text=True)
    (package / "BUILD_REQUIREMENTS.txt").write_text(resolved, encoding="utf-8")
    source_archive(package / "SCATMECH-GUI-source.zip")
    shutil.make_archive(str(DIST / "SCATMECH-GUI-Windows-x64"), "zip", DIST, "SCATMECH-GUI")
    print("\nSUCCESS\nDistribute: " + str(DIST / "SCATMECH-GUI-Windows-x64.zip"))
    print("Recipients extract the complete ZIP and run SCATMECH-GUI\\SCATMECH-GUI.exe.")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as error:
        print(f"\nBUILD FAILED: {error}", file=sys.stderr)
        sys.exit(1)
