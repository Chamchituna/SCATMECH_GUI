from __future__ import annotations

import sys
import traceback

from scatmech_paths import configure_runtime


def main(argv=None):
    args = sys.argv[1:] if argv is None else list(argv)
    for stream in (sys.stdout, sys.stderr):
        if hasattr(stream, "reconfigure"):
            stream.reconfigure(encoding="utf-8")
    if len(args) != 3 or args[0] not in {"--geometry-worker", "--custom-profile-worker"}:
        print("Usage: SCATMECH-Worker --geometry-worker REQUEST OUTPUT\n"
              "   or: SCATMECH-Worker --custom-profile-worker REQUEST OUTPUT_DIR", file=sys.stderr)
        return 2
    try:
        configure_runtime()
        if args[0] == "--geometry-worker":
            from geometry_worker import main as run_worker
        else:
            from custom_profile_worker import main as run_worker
        return run_worker(args[1:])
    except Exception:
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
