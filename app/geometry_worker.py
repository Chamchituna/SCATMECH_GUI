from __future__ import annotations

import json
import os
from pathlib import Path
import sys


def main(argv=None):
    args = sys.argv[1:] if argv is None else argv
    if len(args) != 2:
        print("Expected geometry request and scene output paths.", file=sys.stderr)
        return 2
    if os.name == "posix":
        try:
            os.setsid()
        except OSError:
            if os.getpgrp() != os.getpid():
                print("Could not isolate the geometry worker process.", file=sys.stderr)
                return 2
    try:
        from geometry_scene import build_scene, save_scene
        from geometry_types import GeometryRequest

        payload = json.loads(Path(args[0]).read_text(encoding="utf-8"))
        if not isinstance(payload, dict):
            raise ValueError("Geometry request must be an object.")
        request = GeometryRequest(**payload)
        os.environ["SCATMECH_GEOMETRY_WORKDIR"] = str(Path(args[1]).resolve().parent)

        def progress(message):
            print(json.dumps({"event": "status", "message": str(message)}), flush=True)

        scene = build_scene(request, progress=progress)
        save_scene(scene, args[1])
        progress("Geometry ready.")
        return 0
    except Exception as exc:
        print(f"{type(exc).__name__}: {exc}", file=sys.stderr, flush=True)
        return 1


if __name__ == "__main__":
    sys.exit(main())
