# Third-party components

The GUI's existing MIT license is retained in `LICENSE`; it does not replace
the licenses of bundled third-party components.

- **SCATMECH 7.23:** NIST software and disclaimer in
  `native_src/SCATMECH_LICENSE.txt`. Original source attributions are preserved in
  `native_src/SOURCE_NOTICES.txt`, including the World Scientific
  attribution/permission statement from `tmatrix.cpp`; see also
  `native_src/THIRD_PARTY_NOTICES.txt`. Code comments have been removed.
- **PyQt5:** [GPL v3 or Riverbank commercial license](https://riverbankcomputing.com/software/pyqt/intro).
  This kit installs the PyPI distribution. A proprietary distribution requires
  an appropriate licensing arrangement; this kit does not include one.
- **Qt 5:** supplied through PyQt5-Qt5, with its own
  [Qt licensing terms](https://www.qt.io/licensing/). Retain applicable Qt notices.
- **NumPy, Matplotlib, Python and their dependencies:** each retains its own
  license. The build copies license/notice files found in the installed wheels
  and Python installation into `dist/SCATMECH-GUI/licenses` and records an
  installed dependency inventory. Check that folder for the actual build.
- **PyInstaller:** build tool/bootloader with its
  [GPL exception for generated applications](https://pyinstaller.org/en/stable/license.html).
- **CMake and MSVC:** build tools, not application features; their tool/runtime
  terms apply independently. Native solver builds use the static MSVC runtime.

The output includes `SCATMECH-GUI-source.zip` containing this project's Python
and native build sources, plus `BUILD_REQUIREMENTS.txt` recording the resolved
build environment. Third-party Python/Qt sources are not vendored in this kit.
Review the applicable corresponding-source requirements for your release;
copying notices alone does not fulfill every possible redistribution obligation.
