# SCATMECH GUI

Desktop tools for BRDF, RCW, reflectance, and Mie scattering simulations using
SCATMECH. The GUI includes plotting, geometry previews, and 1D/2D custom-profile
editing. This repository contains the current distribution source, native solver
sources, and the Windows portable-application build kit.

## Run from source on macOS or Linux

Use Python 3.10–3.12, a C++11 compiler, and CMake 3.20 or newer. Python 3.12 is
recommended to match the Windows build environment. Run from the repository root:

```sh
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
cmake -S native_src -B build/native -DCMAKE_BUILD_TYPE=Release
cmake --build build/native --parallel 4
export SCATMECH_BIN="$PWD/build/native/bin"
python app/run_gui.py
```

The native build produces `brdfprog`, `mieprog`, `reflectprog`, `rcwprog`, and
`geometry_export`. The customized drivers supplied here match the GUI's input
and output formats. No separate SCATMECH checkout is required.

## Build the Windows application

Follow [BUILD_WINDOWS.md](BUILD_WINDOWS.md), then run `BUILD_WINDOWS.cmd` on
Windows with CPython 3.12 x64 and Visual Studio Build Tools 2022 installed.
The build compiles the solvers, packages the GUI, and runs an integration check
before creating `dist/SCATMECH-GUI-Windows-x64.zip`.

Recipients extract the complete ZIP and run `SCATMECH-GUI.exe`. This repository
contains source; it does not include prebuilt Windows executables.

## Data and configuration

Source runs save results in `app/DATA/`. The packaged Windows app uses
`%LOCALAPPDATA%\SCATMECH-GUI\DATA`. Generated results, local material files,
plots, caches, virtual environments, and compiled binaries are excluded from
the repository.

| Variable | Purpose |
| --- | --- |
| `SCATMECH_BIN` | Directory containing the five native executables. |
| `SCATMECH_DATA_DIR` | Override the results directory. |
| `SCATMECH_CACHE_DIR` | Override the writable cache directory. |
| `SCATMECH` | Additional search directories for material and grating files. |

## Verify

With the source dependencies installed:

```sh
python -m unittest discover -s tests -v
QT_QPA_PLATFORM=offscreen python app/run_gui.py --smoke-test build/smoke-test.json
```

The integration check requires `SCATMECH_BIN` to point to the compiled solvers.
See [VALIDATION.md](VALIDATION.md) for the checks performed on this snapshot.

## Source and licenses

Application, build, and native code contain no comments; Python docstrings have
also been removed. Runtime strings and C/C++ preprocessor directives are retained.
License and attribution text is preserved in separate notice files.

- GUI: [MIT License](LICENSE), copyright 2026 Woojun Lee.
- Native library: [SCATMECH notice](native_src/SCATMECH_LICENSE.txt),
  [third-party notices](native_src/THIRD_PARTY_NOTICES.txt), and
  [source attributions](native_src/SOURCE_NOTICES.txt).
- Dependencies: [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).
- Provenance and file hashes: [SOURCE_MANIFEST.json](SOURCE_MANIFEST.json).

Woojun Lee · Purdue University · School of Mechanical Engineering
