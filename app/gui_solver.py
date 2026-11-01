import os
import subprocess
from pathlib import Path

from PyQt5.QtCore import QEventLoop, QProcess, QProcessEnvironment, Qt, QTimer
from PyQt5.QtWidgets import QApplication, QProgressDialog


class SolverCancelled(RuntimeError):
    pass


def solver_environment(input, cwd=None):
    environment = os.environ.copy()
    roots = [Path(cwd).resolve()] if cwd is not None else []
    roots.append(Path.cwd())
    search = []
    for line in input.splitlines():
        value = line.partition("!")[0].strip().strip('"')
        if not value:
            continue
        candidate = Path(value).expanduser()
        for root in roots:
            path = candidate if candidate.is_absolute() else root / candidate
            try:
                if path.is_file():
                    search.append(str(path.resolve().parent))
                    break
            except (OSError, ValueError):
                continue
    search.extend(environment.get("SCATMECH", "").split(os.pathsep))
    search = list(dict.fromkeys(path for path in search if path))
    if search:
        environment["SCATMECH"] = os.pathsep.join(search)
    return environment


def run_solver(exe, *, input, parent, cwd=None):

    if getattr(parent, "_solver_running", False):
        raise RuntimeError("A calculation is already running for this module.")
    environment = solver_environment(input, cwd)
    if QApplication.instance() is None:
        return subprocess.run([str(exe)], input=input, cwd=cwd,
                              capture_output=True, text=True, check=False, env=environment,
                              creationflags=subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0)

    process = QProcess()
    process_environment = QProcessEnvironment()
    for name, value in environment.items():
        process_environment.insert(name, value)
    process.setProcessEnvironment(process_environment)
    if cwd is not None:
        process.setWorkingDirectory(str(cwd))
    loop = QEventLoop()
    dialog = QProgressDialog(f"Running {os.path.basename(exe)}…", "Cancel", 0, 0, parent)
    dialog.setWindowTitle("SCATMECH calculation")
    dialog.setWindowModality(Qt.ApplicationModal)
    dialog.setMinimumDuration(0)
    dialog.setAutoClose(False)
    dialog.setAutoReset(False)
    cancelled = False
    failed_to_start = False
    kill_timer = QTimer()
    kill_timer.setSingleShot(True)
    kill_timer.timeout.connect(process.kill)

    def cancel():
        nonlocal cancelled
        cancelled = True
        process.terminate()
        kill_timer.start(1000)

    def on_started():
        if cancelled:
            process.kill()
        else:
            process.write(input.encode("utf-8"))
            process.closeWriteChannel()

    def on_error(error):
        nonlocal failed_to_start
        if error == QProcess.FailedToStart:
            failed_to_start = True
            loop.quit()

    process.started.connect(on_started)
    process.finished.connect(loop.quit)
    process.errorOccurred.connect(on_error)
    dialog.canceled.connect(cancel)
    parent._solver_running = True
    if parent.isVisible():
        dialog.show()
    try:
        process.start(str(exe), [])
        loop.exec_()
        if cancelled:
            raise SolverCancelled("Calculation cancelled.")
        if failed_to_start:
            raise OSError(f"Could not start {exe}: {process.errorString()}")
        stdout = bytes(process.readAllStandardOutput()).decode("utf-8", errors="replace")
        stderr = bytes(process.readAllStandardError()).decode("utf-8", errors="replace")
        returncode = process.exitCode()
        if process.exitStatus() == QProcess.CrashExit and returncode == 0:
            returncode = -1
        return subprocess.CompletedProcess([str(exe)], returncode, stdout, stderr)
    finally:
        kill_timer.stop()
        if process.state() != QProcess.NotRunning:
            process.kill()
            process.waitForFinished(1000)
        dialog.canceled.disconnect(cancel)
        dialog.close()
        dialog.deleteLater()
        parent._solver_running = False
