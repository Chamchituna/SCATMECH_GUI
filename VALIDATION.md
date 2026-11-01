# Distribution validation — November 1, 2026

This source snapshot was checked on macOS with AppleClang 21, CMake 3.31.10,
and Python 3.9.6 using NumPy 2.0.2, Matplotlib 3.9.4, and PyQt5. The documented
Windows build uses CPython 3.12 x64. A Windows executable has not been built or
run for this cleaned snapshot; Windows build and runtime validation remains
required before distributing a Windows binary.

## Cleanup verification

- All 37 Python and PyInstaller-spec files parse with zero comment tokens and
  zero docstrings. Seven Python/spec comments and eleven docstrings were removed.
- Application Python ASTs match the original build-kit application after
  excluding docstrings. The Windows build script additionally includes the
  relocated notices and new documentation in its packaged output.
- All 148 C++ source/header files contain zero comments: 9,247 line comments and
  14 block comments were removed. All 303,523 significant lexer tokens match.
- All 83 C++ translation units produce equivalent preprocessed content,
  ignoring whitespace. Native source line numbers were retained for diagnostics.
- Batch and CMake comments were removed. The Windows manifest remains valid XML.
- Attribution headers and the World Scientific notice are retained in
  `native_src/SOURCE_NOTICES.txt`; native install and Windows packaging include it.

## Build and runtime checks

| Check | Result |
| --- | --- |
| Clean CMake Release build | All five native executables compiled and linked. |
| CMake install | All five executables and all three native notice files installed. |
| Packaging runtime tests | 10 passed. |
| Existing geometry-viewer tests against staged modules | 18 passed. |
| Existing custom-profile-preview tests against staged modules | 6 passed. |
| Existing geometry-worker tests against staged modules | 2 passed. |
| GUI construction | All four pages constructed successfully. |
| BRDF, Reflect, Mie, and RCW | Native calculations, CSV parsing, and Matplotlib plots passed. |
| Geometry worker | Native helper, JSON progress, scene output, and paths containing spaces and Unicode passed. |
| Custom profiles | Both 1D and 2D stdin validation, native geometry, and saved scenes passed. |
| Source archive | All distribution files included; no results, caches, environments, or binaries. |

The native compiler emitted two pre-existing warnings in `axifree.cpp` and
`scateval.cpp`; there were no build errors. The Windows batch launcher and MSVC
manifest behavior were reviewed but not executed on this host.

`SOURCE_MANIFEST.json` contains SHA-256 hashes for every distributed file except
itself. Licenses, documentation, and source are included; local research data,
job configurations, output plots, credentials, and compiled binaries are excluded.
