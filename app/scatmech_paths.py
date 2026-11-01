from __future__ import annotations

import os
import sys
from pathlib import Path
from typing import Final

_DEFAULT_LINUX_DIR: Final[Path] = Path.home() / ".local/share/ME498/SCATMECH"
_DEFAULT_MAC_DIR: Final[Path] = Path.home() / "Local/ME498/SCATMECH"


def is_frozen() -> bool:
    return bool(getattr(sys, "frozen", False))


def get_bundle_dir() -> Path:

    return Path(getattr(sys, "_MEIPASS", Path(__file__).resolve().parent))


def get_scatmech_bin() -> Path:

    override = os.environ.get("SCATMECH_BIN", "").strip()
    if override:
        return Path(override).expanduser().resolve()
    bundled = get_bundle_dir() / "solvers"
    if bundled.is_dir() or is_frozen():
        return bundled
    if sys.platform.startswith("win"):
        local_app_data = os.environ.get("LOCALAPPDATA")
        if local_app_data:
            return Path(local_app_data) / "ME498" / "SCATMECH"
        return _DEFAULT_LINUX_DIR
    if sys.platform == "darwin":
        return _DEFAULT_MAC_DIR
    return _DEFAULT_LINUX_DIR


def get_data_dir() -> Path:
    override = os.environ.get("SCATMECH_DATA_DIR", "").strip()
    if override:
        return Path(override).expanduser().resolve()
    if is_frozen():
        return get_user_dir() / "DATA"
    return Path(__file__).resolve().parent / "DATA"


def get_user_dir() -> Path:

    if sys.platform.startswith("win"):
        base = Path(os.environ.get("LOCALAPPDATA") or Path.home() / "AppData" / "Local")
    elif sys.platform == "darwin":
        base = Path.home() / "Library" / "Application Support"
    else:
        base = Path(os.environ.get("XDG_DATA_HOME") or Path.home() / ".local" / "share")
    return base / "SCATMECH-GUI"


def get_cache_dir() -> Path:
    override = os.environ.get("SCATMECH_CACHE_DIR", "").strip()
    if override:
        return Path(override).expanduser().resolve()
    return get_user_dir() / "cache"


def native_path(path: Path) -> Path:

    path = Path(path).resolve()
    if not sys.platform.startswith("win") or not any(char.isspace() for char in str(path)):
        return path
    import ctypes
    from ctypes import wintypes

    short_path = ctypes.WinDLL("kernel32", use_last_error=True).GetShortPathNameW
    short_path.argtypes = [wintypes.LPCWSTR, wintypes.LPWSTR, wintypes.DWORD]
    short_path.restype = wintypes.DWORD
    size = short_path(str(path), None, 0)
    if size:
        buffer = ctypes.create_unicode_buffer(size)
        length = short_path(str(path), buffer, size)
        if 0 < length < size and not any(char.isspace() for char in buffer.value):
            return Path(buffer.value)
    return path


def get_native_cache_dir() -> Path:

    directory = get_cache_dir() / "custom-profiles"
    if not sys.platform.startswith("win") and not os.environ.get("SCATMECH_CACHE_DIR", "").strip():
        import tempfile
        directory = Path(tempfile.gettempdir()) / ("scatmech-custom-" + str(getattr(os, "getuid", lambda: "user")()))
    directory.mkdir(parents=True, exist_ok=True, mode=0o700)
    directory = native_path(directory)
    if any(char.isspace() for char in str(directory)):
        raise RuntimeError(
            "SCATMECH needs a cache path without spaces and no usable short path is "
            "available. Set SCATMECH_CACHE_DIR to a writable folder without spaces "
            "(for example C:\\SCATMECH-cache on Windows)."
        )
    return directory


def configure_runtime(*, change_directory: bool = False) -> None:
    data = get_data_dir()
    data.mkdir(parents=True, exist_ok=True)
    cache = get_cache_dir()
    cache.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(cache / "matplotlib"))
    search = [str(get_scatmech_bin()), str(data)]
    search.extend(part for part in os.environ.get("SCATMECH", "").split(os.pathsep) if part)
    os.environ["SCATMECH"] = os.pathsep.join(dict.fromkeys(search))
    if change_directory:

        os.chdir(data if is_frozen() else Path(__file__).resolve().parent)


def worker_command(script: Path, *arguments: object) -> tuple[str, list[str]]:

    if is_frozen():
        worker = Path(sys.executable).resolve().with_name("SCATMECH-Worker.exe" if os.name == "nt" else "SCATMECH-Worker")
        flags = {"geometry_worker.py": "--geometry-worker", "custom_profile_worker.py": "--custom-profile-worker"}
        if not worker.is_file():
            raise FileNotFoundError(f"Missing background worker: {worker}. Extract the complete application folder.")
        return str(worker), [flags[Path(script).name], *map(str, arguments)]
    interpreter = Path(sys.executable)
    if interpreter.name.lower() == "pythonw.exe":
        console_interpreter = interpreter.with_name("python.exe")
        if console_interpreter.is_file():
            interpreter = console_interpreter
    return str(interpreter), [str(script), *map(str, arguments)]
