# Native SCATMECH programs

This directory contains all C++ sources required to build the five programs
used by SCATMECH GUI. No SCATMECH checkout, third-party math library, or network
download is needed for this CMake build.

The supported Windows runtime is **64-bit Windows 10 version 1903 or newer,
or Windows 11**. The Windows build requires CMake 3.20 or newer and Visual Studio 2022 Build
Tools with **Desktop development with C++**, including the MSVC x64 toolset
and a Windows SDK. From the repository root:

```bat
cmake -S native_src -B build\native -G "Visual Studio 17 2022" -A x64
cmake --build build\native --config Release --parallel
```

Use a short working path such as `C:\SCATMECH-Build`; the MSVC manifest
input filename must remain within 260 characters after expansion.

The output files are:

```text
build/native/bin/brdfprog.exe
build/native/bin/mieprog.exe
build/native/bin/reflectprog.exe
build/native/bin/rcwprog.exe
build/native/bin/geometry_export.exe
```

The build uses `/MT` for Release to statically link the MSVC C/C++ runtime.
All five Windows programs embed `utf8.manifest`, which selects the UTF-8
process code page so SCATMECH's narrow file APIs and command-line arguments
can use non-ASCII paths passed by the GUI. This setting requires Windows 10
version 1903 or newer. The manifest is merged with the linker's normal
manifest; existing UAC and runtime settings remain in place. See Microsoft's
[UTF-8 process code page documentation](https://learn.microsoft.com/en-us/windows/apps/design/globalizing/use-utf8-code-page).
The programs remain console applications because the GUI communicates with
them using standard input/output. `geometry_export.exe` is prebuilt so users
of the distributed GUI do not need a compiler. The optional CMake install
step copies programs into `bin` and notices into `licenses` beneath the chosen
prefix:

```bat
cmake --install build\native --config Release --prefix build\native-install
```

For a source-level build check on macOS or Linux, use an installed C++11
compiler and omit the Visual Studio generator and architecture flags:

```sh
cmake -S native_src -B build/native -DCMAKE_BUILD_TYPE=Release
cmake --build build/native --parallel 4
```

These produce native host executables, not Windows binaries. Compile and run
the final Windows build on Windows before distribution.

## Source provenance

The upstream library is SCATMECH 7.23, based on the local copy of
[`usnistgov/SCATMECH` commit
`6c5c179ab7bf8de3f5aae64d71c5ad158d7bab4b`](https://github.com/usnistgov/SCATMECH/tree/6c5c179ab7bf8de3f5aae64d71c5ad158d7bab4b),
dated December 18, 2024. Its root-level `.cpp` and `.h` files are included
together with the four command-line drivers. Code comments have been removed
for this distribution; the remaining C++ tokens are unchanged from the build
kit snapshot. Original attribution headers and the World Scientific statement
are preserved in `SOURCE_NOTICES.txt`, with their source-file paths.
Old binaries, object files, backup files, repository metadata, documentation
images, and simulation results are excluded.

The pre-existing local driver adaptations are retained because the GUI
depends on their input/output format:

- `BRDFProg/BRDFProg.cpp` supports the existing optional `--batch` input path.
- `MieProg/MieProg.cpp` accepts theta and phi step sizes and writes the
  corresponding angular grid. Its upstream attribution is retained in
  `SOURCE_NOTICES.txt`.
- `RCWProg/RCWProg.cpp` and `ReflectProg/ReflectProg.cpp` retain their upstream
  code with comments removed.
- `geometry_export.cpp` is copied from this GUI project's `native` directory.

The CMake wrapper is supplied for this distribution kit. Retain
`SCATMECH_LICENSE.txt`, `THIRD_PARTY_NOTICES.txt`, and `SOURCE_NOTICES.txt` with
source and binary distributions. The NIST license includes an explicit
distribution condition for its notice and disclaimer; the third-party notice
originally in `tmatrix.cpp` is reproduced in the separate notice files.
