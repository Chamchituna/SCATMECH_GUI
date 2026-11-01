from __future__ import annotations

import os
import sys
import traceback

from scatmech_paths import configure_runtime, get_user_dir


def main(argv=None):
    args = sys.argv[1:] if argv is None else list(argv)
    if args and args[0] in {"--geometry-worker", "--custom-profile-worker"}:
        from worker_entry import main as worker_main
        return worker_main(args)
    if args and args[0] == "--smoke-test":
        if len(args) != 2:
            return 2
        from package_smoke import run
        return run(args[1])



    for name in ("stdout", "stderr"):
        if getattr(sys, name) is None:
            setattr(sys, name, open(os.devnull, "w", encoding="utf-8"))
    try:
        configure_runtime(change_directory=True)
        from main import QApplication, SCATMECHGui

        app = QApplication([sys.argv[0], *args])
        gui = SCATMECHGui()
        gui.show()
        return app.exec_()
    except Exception:
        details = traceback.format_exc()
        log_message = ""
        try:
            log = get_user_dir() / "startup-error.log"
            log.parent.mkdir(parents=True, exist_ok=True)
            log.write_text(details, encoding="utf-8")
            log_message = f"\n\nDetails saved to {log}"
        except OSError:
            pass
        try:
            from PyQt5.QtWidgets import QApplication, QMessageBox
            app = QApplication.instance() or QApplication([])
            QMessageBox.critical(None, "SCATMECH GUI — Startup Error", details + log_message)
        except Exception:
            sys.stderr.write(details)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
