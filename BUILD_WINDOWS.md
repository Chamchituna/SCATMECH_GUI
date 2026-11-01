# SCATMECH GUI — Windows build kit

Copy this entire folder to a Windows PC, then build a portable application with
`BUILD_WINDOWS.cmd`. This is a source build kit; it does not contain a Windows
executable yet. The source is also available in this repository.

## 1. Install the build prerequisites on Windows

- Windows 10 version 1903 or newer, or Windows 11, **x64**.
- **Python 3.12, 64-bit**, with pip and the Python launcher. The
  [Python 3.12.10 Windows installer](https://www.python.org/downloads/release/python-31210/)
  provides this toolchain. Select **Add python.exe to PATH** and keep pip and
  the Python launcher enabled during installation.
- **Visual Studio Build Tools 2022** (or Visual Studio 2022) with **Desktop
  development with C++**, **MSVC v143 x64/x86 build tools**, and a **Windows SDK**.
  Use the [Microsoft Visual Studio downloads](https://visualstudio.microsoft.com/downloads/)
  or your existing Visual Studio Installer. This script selects the VS 2022
  generator explicitly; a different major Visual Studio release alone is insufficient.
- Internet access for Python package downloads. CMake and PyInstaller are
  installed automatically in this folder's private `.build-venv`.

Recipients of the finished application do not need Python, Visual Studio,
CMake, or a separate SCATMECH installation.

## 2. Build

1. Extract the complete transfer ZIP, preferably into a short local path such as
   `C:\SCATMECH-build`. Do not run it from inside the ZIP viewer.
2. Double-click **BUILD_WINDOWS.cmd**.
3. Wait for **SUCCESS**. The first build compiles the native library and can take
   several minutes. No administrator terminal is needed once prerequisites exist.

The script creates a private Python environment, installs pinned application
dependencies, compiles five native programs, packages the GUI and its worker,
and runs an automated integration check on the actual packaged executable.
The check runs BRDF, Reflect, Mie, RCW, plotting, geometry preview, and both
1D and 2D custom-profile validation, including worker arguments with spaces.

For a terminal without the pause at the end:

```bat
py -3.12 build_windows.py
```

## 3. Distribute

Send **`dist\SCATMECH-GUI-Windows-x64.zip`**. Recipients extract it and launch:

```text
SCATMECH-GUI\SCATMECH-GUI.exe
```

Keep the entire extracted folder together. `SCATMECH-Worker.exe` and `_internal`
are required; sending only `SCATMECH-GUI.exe` will not work. The ZIP also contains
license notices, dependency versions, the build smoke-test report, and the GUI
and native source build kit. This is a portable application, not an installer.

After the build succeeds, open the GUI on Windows and check normal interactive
use, plot navigation, custom-profile editing, and cancelling a long calculation.
Try the final ZIP on a second Windows PC without development tools before
publishing it. macOS validation cannot establish Windows runtime behavior.

## Saved results and inputs

The packaged app saves simulation results under:

```text
%LOCALAPPDATA%\SCATMECH-GUI\DATA
```

Its cache is under `%LOCALAPPDATA%\SCATMECH-GUI\cache`. Startup errors are written
to `%LOCALAPPDATA%\SCATMECH-GUI\startup-error.log`. The app can therefore run from
an installation directory that is not writable. A shortcut may point directly
to the GUI executable; its working directory does not determine the output path.

Optional environment variables:

| Variable | Purpose |
| --- | --- |
| `SCATMECH_DATA_DIR` | Use a different writable results directory. |
| `SCATMECH_CACHE_DIR` | Use a different writable cache directory. |
| `SCATMECH_BIN` | Override the bundled solver directory. Normally leave unset. |
| `SCATMECH` | Additional directories for material/grating input files. |

SCATMECH's native input format has restrictions on whitespace in filenames.
The app uses Windows short paths for generated custom-profile files where
available. If a username contains spaces and the volume has no short aliases,
set `SCATMECH_CACHE_DIR` to a writable path without spaces, such as
`C:\SCATMECH-cache`. User-supplied material/grating files should also use paths
without spaces. These are model inputs; no personal material files, previous
simulation results, plots, or the macOS virtual environment are included here.

## Troubleshooting

- **"Python was not found ... Microsoft Store" / Python 3.12 not found:**
  Windows may be finding its Store shortcut instead of an installed Python.
  Run the [official Python 3.12.10 x64 installer](https://www.python.org/ftp/python/3.12.10/python-3.12.10-amd64.exe),
  select **Add python.exe to PATH**, and keep **pip** and the **Python launcher**
  enabled. Close the failed build window, then reopen `BUILD_WINDOWS.cmd`.
  Confirm the installation in a new Command Prompt with `py -3.12 --version`.
  Disabling a Store alias alone does not install Python. The updated launcher
  checks the actual Python version/architecture and standard installation
  folders even when Python is missing from PATH. No build log is created if
  Python could not start. If you used a custom installation location, run:

  ```bat
  BUILD_WINDOWS.cmd "C:\path\to\Python312\python.exe"
  ```

  This selects the interpreter explicitly. A different Python major/minor
  version or a 32-bit/ARM64 installation will not be used for this x64 kit.
- **Visual Studio 17 2022 not found:** modify/install Build Tools 2022 and its C++
  workload. Restart the command prompt after installation.
- **Build failure:** inspect `build\build.log`; correct the reported issue and
  rerun the batch file. Keep the complete source folder at the same location
  while reusing its build environment.
- **Moved the folder after a build:** remove only the generated `.build-venv`,
  `build`, and `dist` folders, then rebuild. Virtual environments/CMake caches
  store absolute paths. Your app results are stored separately.
- **Smoke test failure:** inspect `build\smoke-test.json` and files in
  `build\smoke-data`. The script does not create a new distributable ZIP when
  its integration check fails. A ZIP from an earlier build may still exist.
- **Missing worker/solver at runtime:** extract the complete distribution again.
  Remove an old `SCATMECH_BIN` override when using the bundled solvers.

## Source and redistribution notices

`app/` is a snapshot of the current GUI with packaging-specific runtime changes.
`native_src/` contains the complete native build inputs, including the customized
BRDF and Mie programs used by this GUI. `SOURCE_MANIFEST.json` records hashes for
the supplied sources. See `native_src/README.md` for upstream provenance.

The existing GUI MIT license is retained in `LICENSE`. Dependency licenses apply
separately. In particular, the PyPI PyQt5 build uses the GPL licensing route;
Riverbank also provides a commercial license. Review
[Riverbank's licensing terms](https://riverbankcomputing.com/software/pyqt/intro)
and `THIRD_PARTY_NOTICES.md` for your intended distribution. The included project
source ZIP is not a substitute for any third-party corresponding-source duties.

See `VALIDATION.md` for checks performed before transferring this kit.
PyInstaller requires building Windows applications on Windows; see its
[platform documentation](https://pyinstaller.org/en/stable/operating-mode.html).
