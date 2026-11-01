from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor, as_completed
from contextlib import contextmanager
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import tempfile

from gui_validation import finite_float, integer
from scatmech_gratings import (
    CROSS_GRATING_SPECS,
    ONE_D_GRATING_SPECS,
    serialize_cross_grating,
    serialize_one_d_grating,
    validate_cross_grating,
    validate_one_d_grating,
)
from scatmech_paths import get_cache_dir, get_scatmech_bin, is_frozen


_ROOT = Path(__file__).resolve().parent
_HELPER = _ROOT / "native" / "geometry_export.cpp"
if not _HELPER.is_file():
    _HELPER = _ROOT.parent / "native_src" / "geometry_export.cpp"
_CACHE = get_cache_dir() / "geometry-preview"
_CROSS_LEAVES = {
    "Generic_CrossGrating", "Cylinder_CrossGrating", "Rectangle_CrossGrating",
    "Sphere_CrossGrating", "Pyramidal_Pit_CrossGrating", "Null_CrossGrating",
}


def _source_dir():
    override = os.environ.get("SCATMECH_SOURCE_DIR") or os.environ.get("SCATMECH_SOURCE")
    source = Path(override).expanduser().resolve() if override else _ROOT.parent / "native_src" / "scatmech" / "code"
    if not (source / "grating.h").is_file() or not (source / "gcross.cpp").is_file():
        raise RuntimeError(f"SCATMECH source was not found at {source}. Set SCATMECH_SOURCE_DIR to the SCATMECH code directory.")
    return source


def _compiler():
    command = shlex.split(os.environ.get("CXX", "c++"))
    executable = shutil.which(command[0]) if command else None
    if not executable:
        raise RuntimeError("A C++ compiler is required for the geometry viewer. Install a C++ toolchain or set CXX.")
    command[0] = str(Path(executable).resolve())
    return command


def _run(command, **kwargs):
    if os.name == "nt":
        kwargs.setdefault("creationflags", subprocess.CREATE_NO_WINDOW)
    result = subprocess.run(command, capture_output=True, text=True, **kwargs)
    if result.returncode:
        detail = (result.stderr or result.stdout).strip()
        raise RuntimeError(f"Geometry helper failed ({result.returncode}):\n{detail[-12000:]}")
    return result


@contextmanager
def _build_lock(path):
    with path.open("a+b") as stream:
        if os.name == "nt":
            import msvcrt
            stream.write(b"0")
            stream.flush()
            stream.seek(0)
            msvcrt.locking(stream.fileno(), msvcrt.LK_LOCK, 1)
        else:
            import fcntl
            fcntl.flock(stream.fileno(), fcntl.LOCK_EX)
        try:
            yield
        finally:
            if os.name == "nt":
                stream.seek(0)
                msvcrt.locking(stream.fileno(), msvcrt.LK_UNLCK, 1)
            else:
                fcntl.flock(stream.fileno(), fcntl.LOCK_UN)


def ensure_exporter(progress=None):
    bundled = get_scatmech_bin() / ("geometry_export.exe" if os.name == "nt" else "geometry_export")
    if bundled.is_file():
        return bundled
    if is_frozen():
        raise RuntimeError(
            f"The packaged geometry helper is missing: {bundled}. "
            "Extract the complete application folder, including _internal/solvers."
        )
    source = _source_dir()
    compiler = _compiler()
    identity = _run(compiler + ["--version"]).stdout
    sources = sorted(source.glob("*.cpp"))
    fingerprint = hashlib.sha256()
    fingerprint.update(str(source).encode())
    fingerprint.update(json.dumps([compiler, identity, "-std=c++11", "-O2"]).encode())
    for path in sorted(source.glob("*.h")) + sources + [_HELPER]:
        fingerprint.update(path.name.encode())
        fingerprint.update(path.read_bytes())
    build = _CACHE / fingerprint.hexdigest()[:24]
    executable = build / ("geometry_export.exe" if os.name == "nt" else "geometry_export")
    if executable.is_file():
        return executable
    build.mkdir(parents=True, exist_ok=True)
    if progress:
        progress("Preparing native SCATMECH geometry helper")
    with _build_lock(build / "build.lock"):
        if executable.is_file():
            return executable
        units = sources + [_HELPER]
        objects = [build / (path.stem + ".o") for path in units]
        pending = [(path, obj) for path, obj in zip(units, objects) if not obj.is_file()]

        def compile_unit(pair):
            path, obj = pair
            temporary = obj.with_suffix(".building.o")
            _run(compiler + ["-std=c++11", "-O2", "-I", str(source), "-c", str(path), "-o", str(temporary)])
            temporary.replace(obj)

        with ThreadPoolExecutor(max_workers=min(4, os.cpu_count() or 2)) as pool:
            tasks = [pool.submit(compile_unit, pair) for pair in pending]
            for done, task in enumerate(as_completed(tasks), start=len(units) - len(pending) + 1):
                task.result()
                if progress:
                    progress(f"Building native geometry helper: {done}/{len(units)} files")
        if progress:
            progress("Linking native geometry helper")
        temporary = executable.with_suffix(".building")
        _run(compiler + [str(path) for path in objects] + ["-o", str(temporary)])
        temporary.replace(executable)
    return executable


def _referenced_files(values):
    roots = [Path.cwd(), get_scatmech_bin()]
    roots.extend(Path(part).expanduser() for part in os.environ.get("SCATMECH", "").split(os.pathsep) if part)
    found = set()
    queue = []
    for value in values:
        path = Path(value.partition("!")[0]).expanduser()
        choices = [path] if path.is_absolute() else [base / path for base in roots]
        for candidate in choices:
            try:
                exists = candidate.is_file()
            except OSError:
                exists = False
            if exists:
                queue.append(candidate.resolve())
                break
    while queue:
        path = queue.pop()
        if path in found:
            continue
        found.add(path)
        try:
            content = path.read_text(encoding="utf-8")
        except (OSError, UnicodeError):
            continue
        for token in re.findall(r"[^\s()\"',;]+", content):
            token = token.partition("!")[0]
            if not token or token in {".", ".."}:
                continue
            filename = Path(token).expanduser()
            choices = [filename] if filename.is_absolute() else [base / filename for base in [path.parent] + roots]
            for candidate in choices:
                try:
                    exists = candidate.is_file()
                except OSError:
                    exists = False
                if exists:
                    resolved = candidate.resolve()
                    if resolved not in found:
                        queue.append(resolved)
                    break
    return sorted(map(str, found))


def referenced_files(node):
    values = []

    def visit(current):
        if not isinstance(current, dict):
            return
        params = current.get("params", {})
        if isinstance(params, dict):
            values.extend(str(value) for value in params.values() if isinstance(value, (str, int, float)))
        children = current.get("children", {})
        if isinstance(children, dict):
            for child in children.values():
                visit(child)

    visit(node)
    return _referenced_files(values)


def export_leaf(node: dict, wavelength_um: float, resolution: int = 48, progress=None) -> dict:
    wavelength = finite_float(wavelength_um, "Preview wavelength", minimum=0, strict_minimum=True)
    resolution = integer(resolution, "Preview resolution", minimum=1, maximum=512)
    model = node.get("model") if isinstance(node, dict) else None
    if model in ONE_D_GRATING_SPECS:
        errors = validate_one_d_grating(node)
        spec = ONE_D_GRATING_SPECS[model]
        serializer = serialize_one_d_grating
    elif model in _CROSS_LEAVES:
        errors = validate_cross_grating(node)
        spec = CROSS_GRATING_SPECS[model]
        serializer = serialize_cross_grating
    else:
        raise ValueError(f"Unsupported native geometry leaf: {model}")
    if errors:
        raise ValueError("\n".join(errors))
    values = serializer(node)[1:]
    fields = [field["name"] for field in spec["fields"] if field["kind"] == "param"]
    files = _referenced_files(values)
    executable = ensure_exporter(progress)
    if progress:
        progress(f"Evaluating {model} at {wavelength:g} µm")
    lines = [model, repr(wavelength), str(resolution)]
    lines.extend(f"{name}\t{value}" for name, value in zip(fields, values))
    environment = os.environ.copy()
    search = [str(Path.cwd()), str(get_scatmech_bin())]
    search.extend(str(Path(path).parent) for path in files)
    search.extend(part for part in environment.get("SCATMECH", "").split(os.pathsep) if part)
    environment["SCATMECH"] = os.pathsep.join(dict.fromkeys(search))
    with tempfile.TemporaryDirectory(prefix="scatmech-geometry-", dir=os.environ.get("SCATMECH_GEOMETRY_WORKDIR")) as folder:
        request = Path(folder) / "request.txt"
        output = Path(folder) / "geometry.json"
        request.write_text("\n".join(lines) + "\n", encoding="utf-8")
        result = _run([str(executable), str(request), str(output)], cwd=folder, env=environment)
        try:
            geometry = json.loads(output.read_text(encoding="utf-8"))
        except (OSError, ValueError) as error:
            raise RuntimeError(f"The geometry helper did not produce valid geometry: {error}") from error
    geometry["files"] = files
    if model in _CROSS_LEAVES and any(
        name in fields and int(values[fields.index(name)]) > resolution for name in ("grid1", "grid2")
    ):
        geometry["notes"].append("The sampled grid preview may omit narrow features; native grid evaluation is unchanged.")
    diagnostics = "\n".join(part.strip() for part in (result.stderr, result.stdout) if part.strip())
    if diagnostics:
        geometry["diagnostics"] = diagnostics[-12000:]
    return geometry
