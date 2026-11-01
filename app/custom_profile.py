from __future__ import annotations

from dataclasses import asdict, dataclass, field, fields
import copy
import hashlib
import json
import math
import os
from pathlib import Path
import re
import tempfile
from typing import Any, Dict, List, Tuple

import numpy as np


@dataclass
class CustomProfileDocument:
    version: int = 1
    dimension: int = 1
    period_um: float = 1.0
    d1_um: float = 0.5
    d2_um: float = 0.5
    zeta_deg: float = 0.0
    film_thickness_um: float = 0.1
    film_material: str = "(1.5,0)"
    medium_i: str = "(1,0)"
    medium_t: str = "(4.05,0.05)"
    nlayers: int = 20
    grid1: int = 1024
    grid2: int = 1024
    curve_segments: int = 192
    points: List[List[float]] = field(default_factory=list)
    features: List[Dict[str, Any]] = field(default_factory=list)
    one_d_mode: str = "surface"
    regions: List[Dict[str, Any]] = field(default_factory=list)
    native_source: str = ""
    native_pstring: str = "(1)"
    native_base_directory: str = ""


@dataclass
class CompiledProfile:
    node: Dict[str, Any]
    faces: np.ndarray
    material_ids: np.ndarray
    materials: List[Tuple[str, str]]
    bounds: np.ndarray
    lattice: np.ndarray
    source_text: str
    fingerprint: str
    files: List[str] = field(default_factory=list)
    one_d_mode: str = "surface"
    cross_section: List[Dict[str, Any]] = field(default_factory=list)
    boundaries_1d: List[Dict[str, Any]] = field(default_factory=list)


def default_profile(dimension, medium_i="(1,0)", medium_t="(4.05,0.05)"):
    if dimension not in (1, 2):
        raise ValueError("Custom profiles support dimensions 1 and 2")
    document = CustomProfileDocument(dimension=dimension, medium_i=medium_i, medium_t=medium_t)
    if dimension == 1:
        document.points = [[-.5, .1], [-.15, .1], [-.1, .3], [.1, .3], [.15, .1], [.5, .1]]
    else:
        document.nlayers = 10
        document.features = [{"id": "feature-1", "name": "Circular hole", "shape": "circle",
                              "role": "hole", "x_um": 0.0, "y_um": 0.0, "rotation_deg": 0.0,
                              "radius_um": .1, "width_um": .2, "height_um": .2,
                              "vertices": [[-.1, -.1], [.1, -.1], [.1, .1], [-.1, .1]],
                              "depth_um": .1, "post_height_um": .1, "material": "(1.5,0)"}]
    return document


def region_profile(preset="core_spacers", medium_i="(1,0)", medium_t="(4.05,0.05)"):
    from custom_profile_1d import region_profile as create
    return create(preset, medium_i, medium_t)


def default_editor_profile(dimension, medium_i="(1,0)", medium_t="(4.05,0.05)"):
    if dimension == 1:
        return region_profile(medium_i=medium_i, medium_t=medium_t)
    return default_profile(dimension, medium_i, medium_t)


def surface_to_regions(document):
    from custom_profile_1d import surface_to_regions as convert
    return convert(document)


def profile_to_dict(document):
    return json.loads(json.dumps(asdict(document), allow_nan=False))


def profile_from_dict(raw):
    if not isinstance(raw, dict):
        raise ValueError("Custom profile must be an object")
    names = {item.name for item in fields(CustomProfileDocument)}
    unknown = set(raw) - names
    if unknown:
        raise ValueError("Unknown custom profile fields: " + ", ".join(sorted(unknown)))
    if raw.get("version", 1) not in (1, 2):
        raise ValueError("Unsupported custom profile version")
    result = CustomProfileDocument(**json.loads(json.dumps(raw)))
    from custom_profile_constraints import ensure_spacer_constraints
    ensure_spacer_constraints(result)
    errors = validate_profile(result)
    if errors:
        raise ValueError("; ".join(errors))
    return result


def profile_lattice(document):
    if document.dimension == 1:
        return np.array([[document.period_um, 0.0], [0.0, document.period_um]], dtype=float)
    angle = math.radians(document.zeta_deg)
    return np.array([[document.d1_um, 0.0],
                     [document.d2_um * math.sin(angle), document.d2_um * math.cos(angle)]], dtype=float)


def _area(vertices):
    return .5 * float(np.sum(vertices[:, 0] * np.roll(vertices[:, 1], -1)
                             - vertices[:, 1] * np.roll(vertices[:, 0], -1)))


def feature_vertices(feature, curve_segments=192):
    shape = feature.get("shape", "circle")
    if shape == "circle":
        angles = np.arange(int(curve_segments)) * (2 * math.pi / int(curve_segments))
        vertices = float(feature.get("radius_um", .1)) * np.column_stack((np.cos(angles), np.sin(angles)))
    elif shape == "rectangle":
        vertices = np.array([[-.5, -.5], [.5, -.5], [.5, .5], [-.5, .5]])
        vertices *= [float(feature.get("width_um", .2)), float(feature.get("height_um", .2))]
    elif shape == "polygon":
        vertices = np.asarray(feature.get("vertices", []), dtype=float)
        if vertices.ndim != 2 or vertices.shape[1:] != (2,) or len(vertices) < 3:
            raise ValueError("Polygons need at least three coordinate pairs")
    else:
        raise ValueError("Unknown feature shape: " + str(shape))
    rotation = math.radians(float(feature.get("rotation_deg", 0.0)))
    matrix = np.array([[math.cos(rotation), math.sin(rotation)], [-math.sin(rotation), math.cos(rotation)]])
    vertices = vertices @ matrix + [float(feature.get("x_um", 0.0)), float(feature.get("y_um", 0.0))]
    return vertices if _area(vertices) >= 0 else vertices[::-1].copy()


def _cross(a, b):
    return float(a[0] * b[1] - a[1] * b[0])


def _point_segment(point, a, b, tolerance):
    edge = b - a
    return (abs(_cross(edge, point - a)) <= tolerance * max(np.linalg.norm(edge), tolerance)
            and np.dot(point - a, point - b) <= tolerance * tolerance)


def _segments_meet(a, b, c, d, tolerance):
    if np.any(np.maximum(np.minimum(a, b), np.minimum(c, d))
              > np.minimum(np.maximum(a, b), np.maximum(c, d)) + tolerance):
        return False
    ab = b - a
    cd = d - c
    values = [_cross(ab, c - a), _cross(ab, d - a), _cross(cd, a - c), _cross(cd, b - c)]
    return (values[0] * values[1] < 0 and values[2] * values[3] < 0
            or any((_point_segment(c, a, b, tolerance), _point_segment(d, a, b, tolerance),
                    _point_segment(a, c, d, tolerance), _point_segment(b, c, d, tolerance))))


def _contains(vertices, point):
    x, y = point
    inside = False
    for a, b in zip(vertices, np.roll(vertices, -1, axis=0)):
        if (a[1] > y) != (b[1] > y):
            crossing = a[0] + (y - a[1]) * (b[0] - a[0]) / (b[1] - a[1])
            if x < crossing:
                inside = not inside
    return inside


def _positive(value):
    try:
        return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value) and value > 0
    except (ValueError, TypeError):
        return False


def _integer(value, minimum, maximum):
    try:
        return not isinstance(value, bool) and int(value) == value and minimum <= value <= maximum
    except (ValueError, TypeError, OverflowError):
        return False


def _material_valid(value):
    return isinstance(value, str) and bool(value.strip()) and not any(char in value for char in "\r\n\x00\"'")


def validate_profile(document):
    errors = []
    if document.version not in (1, 2):
        errors.append("Unsupported custom profile version")
    if document.dimension not in (1, 2):
        return errors + ["Custom profiles support dimensions 1 and 2"]
    if document.one_d_mode not in ("surface", "regions", "native"):
        return errors + ["Unknown 1D editing mode"]
    if document.dimension == 2 and document.one_d_mode != "surface":
        return errors + ["Region and native-source editing require a 1D grating"]
    for key in ("film_material", "medium_i", "medium_t"):
        if not _material_valid(getattr(document, key)):
            errors.append(key + " must contain a single material value or dielectric filename")
    if not _integer(document.nlayers, 1, 10000):
        errors.append("Number of layers must be an integer from 1 to 10000")
    if document.dimension == 1:
        if not _positive(document.period_um):
            return errors + ["Period must be positive and finite"]
        if document.one_d_mode != "surface":
            from custom_profile_1d import validate_one_d
            return errors + validate_one_d(document)
        try:
            points = np.asarray(document.points, dtype=float)
        except (ValueError, TypeError):
            return errors + ["Profile points must be coordinate pairs"]
        if points.ndim != 2 or points.shape[1:] != (2,) or len(points) < 2 or not np.isfinite(points).all():
            return errors + ["Profile needs at least two finite coordinate pairs"]
        tolerance = float(document.period_um) * 1e-10
        if abs(points[0, 0] + document.period_um / 2) > tolerance or abs(points[-1, 0] - document.period_um / 2) > tolerance:
            errors.append("Profile endpoints must lie at the period boundaries")
        if abs(points[0, 1] - points[-1, 1]) > tolerance:
            errors.append("Profile endpoint heights must match for periodic continuity")
        if np.any(np.diff(points[:, 0]) < 0):
            errors.append("Profile points must stay ordered in x; overhangs are unsupported")
        if np.any(points[:, 1] < 0) or not np.any(points[:, 1] > 0):
            errors.append("Profile heights must be nonnegative, with some material above the substrate")
        if np.any(np.linalg.norm(np.diff(points, axis=0), axis=1) <= tolerance):
            errors.append("Adjacent profile points must be distinct")
        for index in range(len(points) - 2):
            a, b, c = points[index:index + 3]
            if abs(a[0] - c[0]) <= tolerance and (b[1] - a[1]) * (c[1] - b[1]) < 0:
                errors.append("Vertical profile steps cannot reverse along the same edge")
                break
        return errors
    if not _positive(document.d1_um) or not _positive(document.d2_um):
        return errors + ["Both lattice periods must be positive and finite"]
    try:
        valid_skew = math.isfinite(float(document.zeta_deg)) and abs(float(document.zeta_deg)) < 89
    except (ValueError, TypeError):
        valid_skew = False
    if not valid_skew:
        return errors + ["Lattice skew must be between -89 and 89 degrees"]
    if not _positive(document.film_thickness_um):
        return errors + ["Film thickness must be positive and finite"]
    for key in ("grid1", "grid2"):
        if not _integer(getattr(document, key), 2, 16384):
            errors.append(key + " must be an integer from 2 to 16384")
    if not _integer(document.curve_segments, 12, 2048):
        return errors + ["Circle segments must be an integer from 12 to 2048"]
    if not isinstance(document.features, list):
        return errors + ["Features must be a list"]
    tolerance = max(float(document.d1_um), float(document.d2_um)) * 1e-10
    inverse = np.linalg.inv(profile_lattice(document))
    polygons = []
    identities = set()
    for index, feature in enumerate(document.features):
        label = "Feature " + str(index + 1)
        if not isinstance(feature, dict):
            errors.append(label + " must be an object")
            continue
        identity = feature.get("id")
        if not isinstance(identity, str) or not identity or identity in identities:
            errors.append(label + " needs a unique ID")
        identities.add(identity if isinstance(identity, str) else str(index))
        if feature.get("role") not in ("hole", "post"):
            errors.append(label + " must be a hole or post")
        if feature.get("shape") == "circle" and not _positive(feature.get("radius_um")):
            errors.append(label + " radius must be positive")
        if feature.get("shape") == "rectangle" and (not _positive(feature.get("width_um")) or not _positive(feature.get("height_um"))):
            errors.append(label + " width and height must be positive")
        if feature.get("role") == "hole":
            depth = feature.get("depth_um")
            if not _positive(depth) or float(depth) > float(document.film_thickness_um) + tolerance:
                errors.append(label + " hole depth must be positive and no greater than the film thickness")
        else:
            if not _positive(feature.get("post_height_um")):
                errors.append(label + " post height must be positive")
            if not _material_valid(feature.get("material", document.film_material)):
                errors.append(label + " post material is invalid")
        try:
            vertices = feature_vertices(feature, document.curve_segments)
        except (ValueError, TypeError, OverflowError) as exc:
            errors.append(label + ": " + str(exc))
            continue
        if not np.isfinite(vertices).all() or abs(_area(vertices)) <= tolerance * tolerance:
            errors.append(label + " must have finite vertices and nonzero area")
            continue
        fractional = vertices @ inverse
        if np.any(np.abs(fractional) >= .5 - 1e-10):
            errors.append(label + " must remain strictly inside the unit cell")
        count = len(vertices)
        edges = list(zip(vertices, np.roll(vertices, -1, axis=0)))
        if any(np.linalg.norm(b - a) <= tolerance for a, b in edges):
            errors.append(label + " has duplicate adjacent vertices")
        intersecting = False
        for first in range(count if feature.get("shape") == "polygon" else 0):
            for second in range(first + 1, count):
                if second == first + 1 or (first == 0 and second == count - 1):
                    continue
                if _segments_meet(*edges[first], *edges[second], tolerance):
                    intersecting = True
                    break
            if intersecting:
                break
        if intersecting:
            errors.append(label + " polygon must not intersect itself")
        for other_index, other in polygons:
            if np.any(np.maximum(np.min(vertices, axis=0), np.min(other, axis=0))
                      > np.minimum(np.max(vertices, axis=0), np.max(other, axis=0)) + tolerance):
                continue
            overlap = (_contains(vertices, other[0]) or _contains(other, vertices[0])
                       or any(_segments_meet(a, b, c, d, tolerance)
                              for a, b in edges for c, d in zip(other, np.roll(other, -1, axis=0))))
            if overlap:
                errors.append(label + " overlaps or touches feature " + str(other_index + 1))
        polygons.append((index, vertices))
    return errors


def _cache_directory():
    from scatmech_paths import get_native_cache_dir
    return get_native_cache_dir()


def _materialize(content, suffix):
    raw = content.encode("utf-8") if isinstance(content, str) else content
    digest = hashlib.sha256(raw).hexdigest()
    directory = _cache_directory()
    target = directory / (digest + suffix)
    if not target.is_file() or target.read_bytes() != raw:
        handle, temporary = tempfile.mkstemp(prefix="profile-", dir=str(directory))
        try:
            with os.fdopen(handle, "wb") as stream:
                stream.write(raw)
            os.replace(temporary, target)
        finally:
            if os.path.exists(temporary):
                os.unlink(temporary)
    return target


def _native_material(value, files):
    value = value.strip()
    numeric = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?"
    if re.fullmatch(r"\(\s*" + numeric + r"\s*,\s*" + numeric + r"\s*\)", value):
        return re.sub(r"\s+", "", value)
    if re.fullmatch(numeric, value):
        return value
    from scatmech_paths import get_scatmech_bin
    name, marker, column = value.partition("!")
    path = Path(name).expanduser()
    candidates = [path] if path.is_absolute() else [Path.cwd() / path, get_scatmech_bin() / path]
    candidates += [Path(root).expanduser() / path for root in os.environ.get("SCATMECH", "").split(os.pathsep) if root]
    for candidate in candidates:
        if candidate.is_file():
            files.append(str(candidate.resolve()))
            safe = _materialize(candidate.read_bytes(), ".nk")
            return str(safe) + (marker + column if marker else "")
    if any(char.isspace() for char in value):
        raise ValueError("Cannot find dielectric file: " + value)
    return value


def native_material_value(value):
    return _native_material(value, [])


def remap_boundary_media(tree, replacements):
    numeric = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?"

    def identity(value):
        if not isinstance(value, str):
            return None
        value = value.strip()
        if re.fullmatch(numeric, value) or re.fullmatch(r"\(\s*" + numeric + r"\s*,\s*" + numeric + r"\s*\)", value):
            return None
        filename, marker, column = value.partition("!")
        return str(Path(filename).expanduser().resolve()) + (marker + column if marker else "")

    changes = {}
    for previous, replacement in replacements:
        for name in ("medium_i", "medium_t"):
            old = previous.get("params", {}).get(name)
            new = replacement.get("params", {}).get(name)
            old_identity, new_identity = identity(old), identity(new)
            if old_identity is None or new_identity is None:
                continue
            if old_identity in changes and identity(changes[old_identity]) != new_identity:
                raise ValueError("Conflicting updated dielectric files for inherited material: " + old)
            changes[old_identity] = new
    result = copy.deepcopy(tree)

    def visit(node):
        if not isinstance(node, dict):
            return
        params = node.get("params", {})
        for name in ("medium_i", "medium_t"):
            key = identity(params.get(name))
            if key in changes:
                params[name] = changes[key]
        for child in node.get("children", {}).values():
            visit(child)

    visit(result)
    return result


def _format(value):
    if abs(float(value)) < 1e-15:
        return "0"
    return format(float(value), ".15g")


def _quad(vertices):
    vertices = np.asarray(vertices, dtype=float)
    if len(vertices) == 3:
        return np.vstack((vertices, vertices[-1]))
    return vertices


def _slab_partition(polygons):
    positions = sorted({-.5, .5, *(round(float(value), 14) for polygon in polygons for value in polygon[:, 0])})
    result = []
    for left, right in zip(positions, positions[1:]):
        if right - left <= 1e-13:
            continue
        middle = (left + right) / 2
        crossings = [(-.5, -.5, -.5), (.5, .5, .5)]
        for polygon in polygons:
            for a, b in zip(polygon, np.roll(polygon, -1, axis=0)):
                if min(a[0], b[0]) < middle < max(a[0], b[0]):
                    slope = (b[1] - a[1]) / (b[0] - a[0])
                    crossings.append((a[1] + (middle - a[0]) * slope,
                                      a[1] + (left - a[0]) * slope,
                                      a[1] + (right - a[0]) * slope))
        crossings.sort(key=lambda item: item[0])
        for bottom, top in zip(crossings, crossings[1:]):
            if top[0] - bottom[0] <= 1e-13:
                continue
            center = (middle, (bottom[0] + top[0]) / 2)
            owner = next((index for index, polygon in enumerate(polygons) if _contains(polygon, center)), None)
            vertices = np.array([[left, bottom[1]], [right, bottom[2]], [right, top[2]], [left, top[1]]])
            for first, second in ((0, 3), (1, 2)):
                if abs(vertices[first, 1] - vertices[second, 1]) <= 1e-12:
                    vertices[first, 1] = vertices[second, 1] = (vertices[first, 1] + vertices[second, 1]) / 2
            result.append((vertices, owner))
    return result


def _native_text(dimension, materials, native_materials, boundaries):
    vertices = []
    lookup = {}
    rows = []
    for polygon, negative, positive in boundaries:
        names = []
        for vertex in polygon:
            key = tuple(_format(component) for component in vertex)
            if key not in lookup:
                lookup[key] = "v" + str(len(vertices))
                vertices.append(lookup[key] + " (" + ",".join(key) + ")")
            names.append(lookup[key])
        if dimension == 2:
            reduced = [name for index, name in enumerate(names) if index == 0 or name != names[index - 1]]
            if len(reduced) > 1 and reduced[0] == reduced[-1]:
                reduced.pop()
            if len(reduced) < 3:
                continue
            names = (["QUAD"] if len(reduced) == 4 else []) + reduced
        rows.append(" ".join(names + [negative, positive]))
    lines = (["CROSSGRATING"] if dimension == 2 else []) + ["PARAMETERS", "profile_version", "END", "MATERIALS"]
    lines += ["m" + str(index) + " " + native_materials[index] for index in range(2, len(materials))]
    lines += ["END", "VERTICES"] + vertices + ["END", "BOUNDARIES"] + rows + ["END", ""]
    return "\n".join(lines)


def compile_profile(document):
    errors = validate_profile(document)
    if errors:
        raise ValueError("; ".join(errors))
    if document.dimension == 1 and document.one_d_mode != "surface":
        from custom_profile_1d import compile_one_d
        return compile_one_d(document)
    lattice = profile_lattice(document)
    materials = [("Incident medium", document.medium_i), ("Substrate", document.medium_t), ("Film", document.film_material)]
    post_ids = {}
    for index, feature in enumerate(document.features if document.dimension == 2 else []):
        if feature["role"] == "post":
            post_ids[index] = len(materials)
            materials.append((feature.get("name", "Post " + str(index + 1)), feature.get("material", document.film_material)))
    files = []
    native_materials = [_native_material(value, files) for label, value in materials]
    names = ["medium_i", "medium_t"] + ["m" + str(index) for index in range(2, len(materials))]
    faces = []
    material_ids = []
    boundaries = []

    def face(vertices, material):
        value = _quad(vertices)
        if np.linalg.norm(np.cross(value[1] - value[0], value[2] - value[0])) + np.linalg.norm(np.cross(value[2] - value[0], value[3] - value[0])) > 1e-18:
            faces.append(value)
            material_ids.append(material)

    if document.dimension == 1:
        points = np.asarray(document.points, dtype=float)
        top = float(np.max(points[:, 1]))
        period = float(document.period_um)
        for a, b in zip(points, points[1:]):
            if a[1] > 0 or b[1] > 0:
                boundaries.append((np.array([a, b]), "medium_i", "m2"))
            if b[0] > a[0]:
                boundaries.append((np.array([[a[0], 0], [b[0], 0]]),
                                   "m2" if a[1] > 0 or b[1] > 0 else "medium_i", "medium_t"))
            face([[a[0], -period / 2, a[1] - top], [b[0], -period / 2, b[1] - top],
                  [b[0], period / 2, b[1] - top], [a[0], period / 2, a[1] - top]], 2 if a[1] > 0 or b[1] > 0 else 1)
            for y in (-period / 2, period / 2):
                face([[a[0], y, -top], [b[0], y, -top], [b[0], y, b[1] - top], [a[0], y, a[1] - top]], 2)
        for point in (points[0], points[-1]):
            face([[point[0], -period / 2, -top], [point[0], period / 2, -top],
                  [point[0], period / 2, point[1] - top], [point[0], -period / 2, point[1] - top]], 2)
        cell = np.array([[-.5, -.5], [.5, -.5], [.5, .5], [-.5, .5]]) @ lattice
    else:
        top = float(document.film_thickness_um) + max((float(feature["post_height_um"]) for feature in document.features if feature["role"] == "post"), default=0.0)
        polygons = [feature_vertices(feature, document.curve_segments) for feature in document.features]
        inverse = np.linalg.inv(lattice)
        partitions = _slab_partition([polygon @ inverse for polygon in polygons])
        for polygon, owner in partitions:
            xy = polygon @ lattice
            feature = document.features[owner] if owner is not None else None
            if feature is None:
                interfaces = [(0.0, 1, 2), (float(document.film_thickness_um), 2, 0)]
            elif feature["role"] == "hole":
                bottom = max(0.0, float(document.film_thickness_um) - float(feature["depth_um"]))
                interfaces = [(0.0, 1, 2), (bottom, 2, 0)] if bottom > 1e-12 else [(0.0, 1, 0)]
            else:
                material = post_ids[owner]
                interfaces = [(0.0, 1, 2), (float(document.film_thickness_um), 2, material),
                              (float(document.film_thickness_um) + float(feature["post_height_um"]), material, 0)]
            for z, below, above in interfaces:
                vertices = np.column_stack((xy, np.full(4, z)))
                boundaries.append((vertices, names[below], names[above]))
                if above == 0:
                    vertices = vertices.copy()
                    vertices[:, 2] -= top
                    face(vertices, below)
        for index, (feature, polygon) in enumerate(zip(document.features, polygons)):
            if feature["role"] == "hole":
                low = max(0.0, float(document.film_thickness_um) - float(feature["depth_um"]))
                high = float(document.film_thickness_um)
                inside, outside = 0, 2
            else:
                low = float(document.film_thickness_um)
                high = low + float(feature["post_height_um"])
                inside, outside = post_ids[index], 0
            for a, b in zip(polygon, np.roll(polygon, -1, axis=0)):
                vertices = np.array([[a[0], a[1], low], [b[0], b[1], low],
                                     [b[0], b[1], high], [a[0], a[1], high]])
                boundaries.append((vertices, names[inside], names[outside]))
                vertices = vertices.copy()
                vertices[:, 2] -= top
                face(vertices, outside if inside == 0 else inside)
        cell = np.array([[-.5, -.5], [.5, -.5], [.5, .5], [-.5, .5]]) @ lattice
        for a, b in zip(cell, np.roll(cell, -1, axis=0)):
            face([[a[0], a[1], -top], [b[0], b[1], -top],
                  [b[0], b[1], document.film_thickness_um - top], [a[0], a[1], document.film_thickness_um - top]], 2)
    substrate_depth = max(float(np.max(np.linalg.norm(lattice, axis=1))) * .08, top * .15)
    face(np.column_stack((cell, np.full(4, -top - substrate_depth))), 1)
    for a, b in zip(cell, np.roll(cell, -1, axis=0)):
        face([[a[0], a[1], -top - substrate_depth], [b[0], b[1], -top - substrate_depth],
              [b[0], b[1], -top], [a[0], a[1], -top]], 1)
    source_text = _native_text(document.dimension, materials, native_materials, boundaries)
    filename = _materialize(source_text, ".grt")
    params = {"medium_i": native_materials[0], "medium_t": native_materials[1],
              "filename": str(filename), "pstring": "(1)", "nlayers": str(int(document.nlayers))}
    if document.dimension == 1:
        model = "Generic_Grating"
        params["period"] = _format(document.period_um)
    else:
        model = "Generic_CrossGrating"
        params.update(d1=_format(document.d1_um), d2=_format(document.d2_um), zeta=_format(document.zeta_deg),
                      grid1=str(int(document.grid1)), grid2=str(int(document.grid2)))
    face_array = np.asarray(faces, dtype=float).reshape(-1, 4, 3)
    bounds = np.column_stack((np.min(face_array, axis=(0, 1)), np.max(face_array, axis=(0, 1))))
    node = {"kind": "grating" if document.dimension == 1 else "cross", "model": model, "params": params}
    if document.dimension == 2:
        node["children"] = {}
    identity = {"document": profile_to_dict(document), "node": node, "source_text": source_text}
    fingerprint = hashlib.sha256(json.dumps(identity, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    return CompiledProfile(node, face_array,
                           np.asarray(material_ids, dtype=int), materials, bounds, lattice,
                           source_text, fingerprint, sorted(set(files)))
