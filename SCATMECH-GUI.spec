from pathlib import Path
import sys

root = Path(SPECPATH)
app = root / "app"
suffix = ".exe" if sys.platform == "win32" else ""
names = ("brdfprog", "mieprog", "reflectprog", "rcwprog", "geometry_export")
binaries = [(str(app / "solvers" / (name + suffix)), "solvers") for name in names]
for source, _ in binaries:
    if not Path(source).is_file():
        raise SystemExit(f"Build the native solvers first: missing {source}")


hidden = sorted(p.stem for p in app.glob("*.py") if p.stem not in {"run_gui", "worker_entry"})
hidden += ["matplotlib.backends.backend_qt5agg", "mpl_toolkits.mplot3d"]
common = dict(pathex=[str(app)], binaries=binaries, datas=[],
              hiddenimports=hidden, hookspath=[], runtime_hooks=[],
              excludes=["tkinter", "PyQt6", "PySide2", "PySide6"],
              hooksconfig={"matplotlib": {"backends": ["Qt5Agg", "Agg"]}}, noarchive=False)
gui = Analysis([str(app / "run_gui.py")], **common)
worker = Analysis([str(app / "worker_entry.py")], **common)
gui_exe = EXE(PYZ(gui.pure), gui.scripts, [], exclude_binaries=True,
              name="SCATMECH-GUI", console=False, debug=False,
              strip=False, upx=False, contents_directory="_internal")
worker_exe = EXE(PYZ(worker.pure), worker.scripts, [], exclude_binaries=True,
                 name="SCATMECH-Worker", console=True, debug=False,
                 strip=False, upx=False, contents_directory="_internal")
bundle = COLLECT(gui_exe, worker_exe, gui.binaries, gui.datas,
                 worker.binaries, worker.datas, name="SCATMECH-GUI",
                 strip=False, upx=False)
