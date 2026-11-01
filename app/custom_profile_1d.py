from __future__ import annotations

import copy
import hashlib
import json
import math
import os
from pathlib import Path
import re

import numpy as np


def region_profile(preset="core_spacers", medium_i="(1,0)", medium_t="(4.05,0.05)"):
    from custom_profile import CustomProfileDocument
    document = CustomProfileDocument(version=2, one_d_mode="regions", medium_i=medium_i, medium_t=medium_t)
    film = {"id": "film", "name": "Film", "material": "(1.46,0)",
            "vertices": [[-.5, 0], [.5, 0], [.5, .05], [-.5, .05]]}
    core = {"id": "core", "name": "Core", "material": "(3.5,0.02)",
            "vertices": [[-.18, .05], [.18, .05], [.14, .1], [.11, .45],
                         [.04, .5], [-.04, .5], [-.11, .45], [-.14, .1]]}
    if preset == "core_spacers":
        left = {"id": "left-spacer", "name": "Left spacer", "material": "(1.8,0)",
                "vertices": [[-.38, .05], [-.18, .05], [-.14, .1],
                             [-.12285714285714286, .3], [-.31, .2]]}
        right = {"id": "right-spacer", "name": "Right spacer", "material": "(1.8,0)",
                 "vertices": [[-x, z] for x, z in reversed(left["vertices"])]}
        document.regions = [film, left, right, core]
        from custom_profile_constraints import ensure_spacer_constraints
        ensure_spacer_constraints(document)
    elif preset in ("overhang", "void"):
        core["vertices"] = [[-.1, .05], [.1, .05], [.1, .25], [.3, .25],
                            [.3, .4], [-.3, .4], [-.3, .25], [-.1, .25]]
        document.regions = [film, core]
        if preset == "void":
            document.regions.append({"id": "void", "name": "Enclosed air cavity", "material": "@incident",
                                     "vertices": [[-.06, .11], [.06, .11], [.06, .2], [-.06, .2]]})
    elif preset == "empty":
        document.regions = []
    else:
        raise ValueError("Unknown 1D region preset: " + str(preset))
    return document


def surface_to_regions(document):
    from custom_profile import validate_profile
    result = copy.deepcopy(document)
    if document.dimension != 1:
        raise ValueError("Only 1D profiles can be converted to regions")
    if document.one_d_mode == "regions":
        return result
    if document.one_d_mode != "surface":
        raise ValueError("Native boundary source cannot be converted to editable regions")
    errors = validate_profile(document)
    if errors:
        raise ValueError("; ".join(errors))
    points = np.asarray(document.points, dtype=float)
    result.version = 2
    result.one_d_mode = "regions"
    result.regions = []
    pieces = []
    current = []
    for point in points.tolist():
        if not current or point != current[-1]:
            current.append(point)
        if point[1] == 0 and len(current) > 1:
            pieces.append(current)
            current = [point]
    if len(current) > 1:
        pieces.append(current)
    for piece in pieces:
        if max(point[1] for point in piece) <= 0:
            continue
        vertices = [[piece[0][0], 0], [piece[-1][0], 0]] + list(reversed(piece))
        vertices = _clean_polygon(vertices, document.period_um * 1e-12)
        if len(vertices) >= 3:
            index = len(result.regions) + 1
            result.regions.append({"id": "surface-" + str(index), "name": "Surface profile " + str(index),
                                   "material": document.film_material, "vertices": vertices.tolist()})
    return result


def _clean_polygon(vertices, tolerance):
    result = []
    for point in np.asarray(vertices, dtype=float):
        if not result or np.linalg.norm(point - result[-1]) > tolerance:
            result.append(point)
    if len(result) > 1 and np.linalg.norm(result[-1] - result[0]) <= tolerance:
        result.pop()
    changed = True
    while changed and len(result) > 3:
        changed = False
        for index in range(len(result)):
            a, b, c = result[index - 1], result[index], result[(index + 1) % len(result)]
            cross = np.cross(b - a, c - b)
            if abs(cross) <= tolerance * max(np.linalg.norm(c - a), tolerance) and np.dot(b - a, c - b) >= 0:
                result.pop(index)
                changed = True
                break
    return np.asarray(result, dtype=float).reshape(-1, 2)


def validate_one_d(document):
    from custom_profile import _area, _material_valid, _segments_meet
    errors = []
    if document.one_d_mode == "native":
        if not isinstance(document.native_source, str) or not document.native_source.strip():
            return ["Enter a complete Generic_Grating boundary description"]
        if "\x00" in document.native_source:
            errors.append("Native source cannot contain NUL characters")
        content = re.sub(r";[^\n]*", "", document.native_source)
        for name in ("PARAMETERS", "MATERIALS", "VERTICES", "BOUNDARIES"):
            if re.search(r"\b" + name + r"\b", content) is None:
                errors.append("Native source is missing the " + name + " block")
        if not isinstance(document.native_pstring, str) or not document.native_pstring.strip() or any(char in document.native_pstring for char in "\r\n\x00"):
            errors.append("Native parameter values must be a nonempty single-line expression")
        if not isinstance(document.native_base_directory, str) or "\x00" in document.native_base_directory:
            errors.append("Native source directory must be a valid path")
        return errors
    if not isinstance(document.regions, list):
        return ["Regions must be a list"]
    identities = set()
    tolerance = float(document.period_um) * 1e-11
    for index, region in enumerate(document.regions):
        label = "Region " + str(index + 1)
        if not isinstance(region, dict):
            errors.append(label + " must be an object")
            continue
        identity = region.get("id")
        if not isinstance(identity, str) or not identity or identity in identities:
            errors.append(label + " needs a unique ID")
        if isinstance(identity, str):
            identities.add(identity)
        material = region.get("material")
        if material not in ("@incident", "@substrate") and not _material_valid(material):
            errors.append(label + " needs a material value or dielectric filename")
        try:
            vertices = np.asarray(region.get("vertices", []), dtype=float)
        except (ValueError, TypeError):
            errors.append(label + " vertices must be coordinate pairs")
            continue
        if vertices.ndim != 2 or vertices.shape[1:] != (2,) or len(vertices) < 3 or not np.isfinite(vertices).all():
            errors.append(label + " needs at least three finite x–z vertices")
            continue
        if abs(_area(vertices)) <= tolerance * tolerance:
            errors.append(label + " must enclose a nonzero area")
        edges = list(zip(vertices, np.roll(vertices, -1, axis=0)))
        if any(np.linalg.norm(b - a) <= tolerance for a, b in edges):
            errors.append(label + " has duplicate adjacent vertices")
        if any(_segments_meet(*edges[first], *edges[second], tolerance)
               for first in range(len(edges)) for second in range(first + 1, len(edges))
               if second != first + 1 and (first != 0 or second != len(edges) - 1)):
            errors.append(label + " polygon must not intersect itself")
        for index in range(len(vertices)):
            a, b, c = vertices[index - 1], vertices[index], vertices[(index + 1) % len(vertices)]
            if abs(np.cross(b - a, c - b)) <= tolerance * max(np.linalg.norm(c - a), tolerance) and np.dot(b - a, c - b) < 0:
                errors.append(label + " has a reversed edge")
                break
    if not errors:
        from custom_profile_constraints import constraint_errors
        errors.extend(constraint_errors(document))
    return errors


def _cross(a, b):
    return float(a[0] * b[1] - a[1] * b[0])


def _segment_cuts(a, b, c, d, tolerance):
    delta, other = b - a, d - c
    denominator = _cross(delta, other)
    if abs(denominator) > tolerance * max(np.linalg.norm(delta), np.linalg.norm(other), tolerance):
        first = _cross(c - a, other) / denominator
        second = _cross(c - a, delta) / denominator
        if -tolerance <= first <= 1 + tolerance and -tolerance <= second <= 1 + tolerance:
            return [min(1., max(0., first))]
        return []
    length = np.dot(delta, delta)
    if abs(_cross(c - a, delta)) <= tolerance * max(math.sqrt(length), tolerance):
        return [min(1., max(0., value)) for value in (np.dot(c - a, delta) / length, np.dot(d - a, delta) / length)
                if -tolerance <= value <= 1 + tolerance]
    return []


def _clip_segment(a, b, left, right, tolerance):
    delta = b - a
    if abs(delta[0]) <= tolerance:
        return (a.copy(), b.copy()) if left - tolerance <= a[0] <= right + tolerance else None
    values = sorted(((left - a[0]) / delta[0], (right - a[0]) / delta[0]))
    low, high = max(0., values[0]), min(1., values[1])
    if high - low <= tolerance:
        return None
    return a + low * delta, a + high * delta


def _distance(point, a, b):
    edge = b - a
    factor = np.clip(np.dot(point - a, edge) / np.dot(edge, edge), 0., 1.)
    return float(np.linalg.norm(point - a - factor * edge))


def _arrangement(document, material_ids):
    from custom_profile import _contains
    period = float(document.period_um)
    left, right = -period / 2, period / 2
    tolerance = period * 1e-11
    copies = []
    segments = [(np.array([left, 0.]), np.array([right, 0.]))]
    for region, material in zip(document.regions, material_ids):
        polygon = np.asarray(region["vertices"], dtype=float)
        low = math.ceil((left - float(polygon[:, 0].max())) / period)
        high = math.floor((right - float(polygon[:, 0].min())) / period)
        if high - low > 4096:
            raise ValueError("A region spans too many periods to evaluate; reduce its horizontal extent")
        for shift in range(low, high + 1):
            translated = polygon + [shift * period, 0.]
            copies.append((translated, material))
            for a, b in zip(translated, np.roll(translated, -1, axis=0)):
                clipped = _clip_segment(a, b, left, right, tolerance)
                if clipped is not None and np.linalg.norm(clipped[1] - clipped[0]) > tolerance:
                    segments.append(clipped)
    if len(segments) > 6000:
        raise ValueError("The periodic design has too many boundary segments to evaluate interactively")

    def classify(point):
        value = np.asarray(point, dtype=float).copy()
        value[0] = (value[0] - left) % period + left
        for polygon, material in reversed(copies):
            if _contains(polygon, value):
                return material
        return 1 if value[1] < 0 else 0

    split = []
    for index, (a, b) in enumerate(segments):
        cuts = [0., 1.]
        for other, (c, d) in enumerate(segments):
            if other != index:
                cuts.extend(_segment_cuts(a, b, c, d, tolerance))
        cuts = sorted(set(round(value, 13) for value in cuts))
        for start, stop in zip(cuts, cuts[1:]):
            first, second = a + start * (b - a), a + stop * (b - a)
            if np.linalg.norm(second - first) > tolerance:
                split.append((first, second))
    seen = set()
    interfaces = []
    coordinates = []
    for axis in range(2):
        values = sorted({float(point[axis]) for segment in split for point in segment})
        preferred = [float(point[axis]) for segment in segments for point in segment]
        groups = []
        for value in values:
            if not groups or value - groups[-1][0] > tolerance:
                groups.append([value])
            else:
                groups[-1].append(value)
        mapping = {}
        for group in groups:
            representative = next((value for value in preferred if group[0] - tolerance <= value <= group[-1] + tolerance), group[0])
            mapping.update((value, representative) for value in group)
        coordinates.append(mapping)

    def weld(point):
        return np.array([coordinates[axis][float(point[axis])] for axis in range(2)])

    for a, b in split:
        a, b = weld(a), weld(b)
        if np.linalg.norm(b - a) <= tolerance:
            continue
        if abs(a[0] - b[0]) <= tolerance and abs(a[0] - right) <= tolerance:
            continue
        key = tuple(sorted((tuple(a), tuple(b))))
        if key in seen:
            continue
        seen.add(key)
        delta = b - a
        length = np.linalg.norm(delta)
        normal = np.array([-delta[1], delta[0]]) / length
        midpoint = (a + b) / 2
        distances = [_distance(midpoint, c, d) for c, d in segments]
        clearance = min((distance for distance in distances if distance > tolerance), default=period)
        offset = min(length * 1e-5, period * 1e-7, clearance * .1)
        first, second = classify(midpoint + normal * offset), classify(midpoint - normal * offset)
        if first != second:
            interfaces.append((np.array([a, b]), first, second))
    if not interfaces:
        interfaces = [(np.array([[left, 0.], [right, 0.]]), 0, 1)]
    return interfaces, classify


def _partition(interfaces, classify, left, right, bottom, top, tolerance):
    positions = [left, right]
    for segment, first, second in interfaces:
        positions.extend(segment[:, 0])
    positions.sort()
    positions = [value for index, value in enumerate(positions) if index == 0 or value - positions[index - 1] > tolerance]
    cells = []
    for start, stop in zip(positions, positions[1:]):
        middle = (start + stop) / 2
        crossings = [(bottom, bottom, bottom), (top, top, top)]
        for segment, first, second in interfaces:
            a, b = segment
            if min(a[0], b[0]) < middle < max(a[0], b[0]):
                slope = (b[1] - a[1]) / (b[0] - a[0])
                crossings.append((a[1] + (middle - a[0]) * slope,
                                  a[1] + (start - a[0]) * slope,
                                  a[1] + (stop - a[0]) * slope))
        crossings.sort(key=lambda value: value[0])
        for low, high in zip(crossings, crossings[1:]):
            if high[0] - low[0] <= tolerance:
                continue
            polygon = _clean_polygon([[start, low[1]], [stop, low[2]], [stop, high[2]], [start, high[1]]], tolerance)
            if len(polygon) >= 3:
                cells.append({"vertices": (polygon - [0., top]).tolist(),
                              "material_id": int(classify([middle, (low[0] + high[0]) / 2]))})
    return cells


def compile_one_d(document):
    from custom_profile import CompiledProfile, _format, _materialize, _native_material, _native_text, profile_lattice, profile_to_dict
    materials = [("Incident medium", document.medium_i), ("Substrate", document.medium_t)]
    files = []
    native = [_native_material(value, files) for label, value in materials]
    lattice = profile_lattice(document)
    faces, face_ids, cells, boundary_data = [], [], [], []
    period = float(document.period_um)
    left, right = -period / 2, period / 2
    if document.one_d_mode == "native":
        source = _native_source(document, files)
        pstring = document.native_pstring
        bounds = np.array([[left, right], [left, right], [-period * .2, 0.]])
    else:
        region_ids = []
        for index, region in enumerate(document.regions):
            value = region["material"]
            if value in ("@incident", "@substrate"):
                region_ids.append(0 if value == "@incident" else 1)
            else:
                region_ids.append(len(materials))
                materials.append((region.get("name") or "Region " + str(index + 1), value))
                native.append(_native_material(value, files))
        interfaces, classify = _arrangement(document, region_ids)
        top = max(float(segment[:, 1].max()) for segment, first, second in interfaces)
        low = min(float(segment[:, 1].min()) for segment, first, second in interfaces)
        bottom = low - max(.08 * period, .15 * (top - low))
        cells = _partition(interfaces, classify, left, right, bottom, top, period * 1e-11)
        names = ["medium_i", "medium_t"] + ["m" + str(index) for index in range(2, len(materials))]
        source = _native_text(1, materials, native, [(segment, names[first], names[second]) for segment, first, second in interfaces])
        pstring = "(1)"
        for cell in cells:
            vertices = np.asarray(cell["vertices"])
            if len(vertices) == 3:
                vertices = np.vstack((vertices, vertices[-1]))
            for y in (left, right):
                faces.append(np.column_stack((vertices[:, 0], np.full(4, y), vertices[:, 1])))
                face_ids.append(cell["material_id"])
        for segment, first, second in interfaces:
            a, b = segment - [0., top]
            boundary_data.append({"start": a.tolist(), "end": b.tolist(), "left": first, "right": second})
            faces.append([[a[0], left, a[1]], [b[0], left, b[1]], [b[0], right, b[1]], [a[0], right, a[1]]])
            face_ids.append(first if first != 0 else second)
        bounds = np.array([[left, right], [left, right], [bottom - top, 0.]])
    filename = _materialize(source, ".grt")
    node = {"kind": "grating", "model": "Generic_Grating", "params": {
        "period": _format(period), "medium_i": native[0], "medium_t": native[1], "filename": str(filename),
        "pstring": pstring, "nlayers": str(int(document.nlayers))}}
    identity = {"document": profile_to_dict(document), "node": node, "source_text": source}
    fingerprint = hashlib.sha256(json.dumps(identity, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    return CompiledProfile(node, np.asarray(faces, dtype=float).reshape(-1, 4, 3), np.asarray(face_ids, dtype=int),
                           materials, bounds, lattice, source, fingerprint, sorted(set(files)), document.one_d_mode, cells, boundary_data)


def _native_source(document, files):
    from custom_profile import _materialize
    from scatmech_paths import get_scatmech_bin
    base = Path(document.native_base_directory).expanduser() if document.native_base_directory else Path.cwd()
    roots = [base, Path.cwd(), get_scatmech_bin()]
    roots.extend(Path(value).expanduser() for value in os.environ.get("SCATMECH", "").split(os.pathsep) if value)
    cached = {}
    active = set()

    def rewrite_tables(text, search):
        def table(match):
            replacement = resolve(match.group(1), search)
            return "@" + (Path(replacement).name if replacement else match.group(1)) + "("

        result = []
        for line in text.splitlines(keepends=True):
            body, separator, comment = line.partition(";")
            result.append(re.sub(r"@([^\s()]+)\(", table, body) + separator + comment)
        return "".join(result)

    def cache(candidate):
        identity = candidate.resolve()
        if identity in active:
            raise ValueError("Cyclic native table dependency: " + str(identity))
        if identity in cached:
            return cached[identity]
        files.append(str(identity))
        active.add(identity)
        try:
            content = candidate.read_bytes()
            try:
                text = content.decode("utf-8")
            except UnicodeDecodeError:
                text = None
            if text is not None:
                content = rewrite_tables(text, [candidate.parent] + roots).encode("utf-8")
            cached[identity] = _materialize(content, ".nk")
            return cached[identity]
        finally:
            active.remove(identity)

    def resolve(value, search=None):
        filename, marker, column = value.partition("!")
        path = Path(filename).expanduser()
        for candidate in ([path] if path.is_absolute() else [root / path for root in (search if search is not None else roots)]):
            try:
                exists = candidate.is_file()
            except OSError:
                exists = False
            if exists:
                return str(cache(candidate)) + (marker + column if marker else "")
        return None

    source = rewrite_tables(document.native_source, roots)
    tokens = []
    start = None
    depth = 0
    quote = ""
    index = 0
    while index < len(source):
        char = source[index]
        if start is None:
            if char == ";":
                stop = source.find("\n", index)
                index = len(source) if stop < 0 else stop
                continue
            if char.isspace():
                index += 1
                continue
            start = index
        if char in "\"'":
            if not quote:
                quote = char
            elif quote == char:
                quote = ""
        elif not quote:
            if char == "(":
                depth += 1
            elif char == ")":
                depth -= 1
            if depth == 0 and (char.isspace() or char == ";"):
                tokens.append((start, index, source[start:index]))
                start = None
                continue
        index += 1
    if start is not None:
        tokens.append((start, len(source), source[start:]))
    material_start = next((index + 1 for index, token in enumerate(tokens) if token[2] == "MATERIALS"), None)
    if material_start is None:
        return source
    replacements = []
    index = material_start
    while index < len(tokens) and tokens[index][2] != "END":
        name = tokens[index][2]
        count = {"ANISO": 3, "MAGNETIC": 2, "ANISOMAGNETIC": 6}.get(name, 1)
        index += 2 if name in ("ANISO", "MAGNETIC", "ANISOMAGNETIC") else 1
        for start, stop, value in tokens[index:index + count]:
            literal = value[1:-1] if len(value) > 1 and value[0] == value[-1] and value[0] in "\"'" else value
            if "(" not in literal and ")" not in literal:
                replacement = resolve(literal)
                if replacement:
                    replacements.append((start, stop, replacement))
        index += count
    for start, stop, replacement in reversed(replacements):
        source = source[:start] + replacement + source[stop:]
    return source
