import math

import numpy as np


def _clip(polygon, axis, limit, keep_above):
    if not len(polygon):
        return polygon
    result = []
    previous = polygon[-1]
    previous_inside = previous[axis] >= limit - 1e-12 if keep_above else previous[axis] <= limit + 1e-12
    for current in polygon:
        current_inside = current[axis] >= limit - 1e-12 if keep_above else current[axis] <= limit + 1e-12
        if current_inside != previous_inside:
            fraction = (limit - previous[axis]) / (current[axis] - previous[axis])
            result.append(previous + fraction * (current - previous))
        if current_inside:
            result.append(current)
        previous, previous_inside = current, current_inside
    return np.asarray(result, dtype=float).reshape(-1, 3)


def _area(polygon):
    return sum(np.linalg.norm(np.cross(polygon[index] - polygon[0], polygon[index + 1] - polygon[0])) / 2 for index in range(1, len(polygon) - 1))


def _quads(polygon):
    if len(polygon) == 4:
        return [polygon]
    return [np.array([polygon[0], polygon[index], polygon[index + 1], polygon[index + 1]]) for index in range(1, len(polygon) - 1)]


def _subdivide(face, spacing):
    nu = max(1, int(math.ceil(max(np.linalg.norm(face[1] - face[0]), np.linalg.norm(face[2] - face[3])) / spacing)))
    nv = max(1, int(math.ceil(max(np.linalg.norm(face[3] - face[0]), np.linalg.norm(face[2] - face[1])) / spacing)))
    u = np.linspace(0, 1, min(nu, 96) + 1)
    v = np.linspace(0, 1, min(nv, 96) + 1)
    a = (1 - u[:, None]) * face[0] + u[:, None] * face[1]
    b = (1 - u[:, None]) * face[3] + u[:, None] * face[2]
    points = (1 - v[None, :, None]) * a[:, None] + v[None, :, None] * b[:, None]
    return np.stack((points[:-1, :-1], points[1:, :-1], points[1:, 1:], points[:-1, 1:]), axis=2).reshape(-1, 4, 3)


def _plane_segment(polygon, axis, position):
    points = []
    for a, b in zip(polygon, np.roll(polygon, -1, axis=0)):
        da, db = a[axis] - position, b[axis] - position
        if abs(da) <= 1e-10:
            points.append(a)
        if da * db < -1e-20:
            points.append(a + (-da / (db - da)) * (b - a))
    if len(points) < 2:
        return None
    other = 1 - axis
    points = np.asarray(points)
    a, b = points[np.argmin(points[:, other])], points[np.argmax(points[:, other])]
    if b[other] - a[other] <= 1e-12:
        return None
    return a[[other, 2]], b[[other, 2]]


def _caps(polygons, ids, limits, base, bottom, film_level):
    faces, materials = [], []
    for axis in (0, 1):
        other = 1 - axis
        for position in limits[axis]:
            segments = []
            for polygon, material in zip(polygons, ids):
                if np.all(abs(polygon[:, axis] - position) < 1e-10):
                    continue
                segment = _plane_segment(polygon, axis, position)
                if segment is not None:
                    a, b = segment
                    segments.append((a, b, material))
            positions = sorted({round(float(value), 12) for a, b, material in segments for value in (a[0], b[0])} | set(limits[other]))
            for left, right in zip(positions, positions[1:]):
                if right - left <= 1e-12:
                    continue
                middle = (left + right) / 2
                if not limits[other][0] <= middle <= limits[other][1]:
                    continue
                candidates = []
                for a, b, material in segments:
                    if a[0] - 1e-10 <= middle <= b[0] + 1e-10:
                        slope = (b[1] - a[1]) / (b[0] - a[0])
                        candidates.append((a[1] + slope * (middle - a[0]), a[1] + slope * (left - a[0]), a[1] + slope * (right - a[0]), material))
                if not candidates:
                    continue
                _, zleft, zright, material = max(candidates, key=lambda row: row[0])
                spans = [(bottom, bottom, base, base, 1)]
                if material == 2:
                    spans.append((base, base, zleft, zright, 2))
                elif material >= 3:
                    spans.append((base, base, film_level, film_level, 2))
                    spans.append((film_level, film_level, zleft, zright, material))
                for low_left, low_right, high_left, high_right, surface_material in spans:
                    if max(high_left - low_left, high_right - low_right) <= 1e-12:
                        continue
                    face = np.zeros((4, 3))
                    face[:, axis] = position
                    face[:, other] = [left, right, right, left]
                    face[:, 2] = [low_left, low_right, high_right, high_left]
                    faces.append(face)
                    materials.append(surface_material)
    return faces, materials


def _region_design_mesh(compiled, repeats):
    nx, ny = (max(1, min(5, int(value))) for value in repeats)
    period = float(compiled.lattice[0, 0])
    left, right = -nx * period / 2, nx * period / 2
    front, back = -ny * period / 2, ny * period / 2
    cells = compiled.cross_section
    boundaries = compiled.boundaries_1d
    faces, identifiers = [], []
    tolerance = period * 1e-10

    def emit(polygon, material):
        if len(polygon) < 3 or _area(polygon) <= tolerance * tolerance:
            return
        for face in _quads(np.asarray(polygon, dtype=float)):
            pieces = _subdivide(face, period / 12)
            faces.extend(pieces)
            identifiers.extend([material] * len(pieces))

    tiles = range(math.ceil(-nx / 2 - .5 + 1e-10), math.floor(nx / 2 + .5 - 1e-10) + 1)
    for tile in tiles:
        shift = tile * period
        for boundary in boundaries:
            negative, positive = boundary["left"], boundary["right"]
            if negative != 0 and positive != 0:
                continue
            material = positive if negative == 0 else negative
            if material == 0:
                continue
            a, b = np.asarray(boundary["start"]) + [shift, 0], np.asarray(boundary["end"]) + [shift, 0]
            if any(abs(a[0] - edge) <= tolerance and abs(b[0] - edge) <= tolerance for edge in (left, right)):
                continue
            polygon = np.array([[a[0], front, a[1]], [b[0], front, b[1]],
                                [b[0], back, b[1]], [a[0], back, a[1]]])
            polygon = _clip(_clip(polygon, 0, left, True), 0, right, False)
            emit(polygon, material)
        for cell in cells:
            material = cell["material_id"]
            if not material:
                continue
            section = np.asarray(cell["vertices"], dtype=float) + [shift, 0]
            polygon = np.column_stack((section[:, 0], np.zeros(len(section)), section[:, 1]))
            polygon = _clip(_clip(polygon, 0, left, True), 0, right, False)
            if len(polygon) < 3 or _area(polygon) <= tolerance * tolerance:
                continue
            for y in (front, back):
                cap = polygon.copy()
                cap[:, 1] = y
                emit(cap, material)
            for x in (left, right):
                points = polygon[np.abs(polygon[:, 0] - x) <= tolerance, 2]
                if len(points) >= 2 and np.ptp(points) > tolerance:
                    low, high = float(np.min(points)), float(np.max(points))
                    emit(np.array([[x, front, low], [x, back, low], [x, back, high], [x, front, high]]), material)
    bottom = min(point[1] for cell in cells for point in cell["vertices"])
    emit(np.array([[left, front, bottom], [right, front, bottom], [right, back, bottom], [left, back, bottom]]), 1)
    faces = np.asarray(faces, dtype=float).reshape(-1, 4, 3)
    bounds = np.column_stack((faces.min(axis=(0, 1)), faces.max(axis=(0, 1))))
    return faces, np.asarray(identifiers, dtype=int), compiled.materials, bounds


def build_design_mesh(compiled, repeats=(3, 3)):
    if getattr(compiled, "one_d_mode", "surface") == "regions":
        return _region_design_mesh(compiled, repeats)
    nx, ny = (max(1, min(5, int(value))) for value in repeats)
    lattice = np.asarray(compiled.lattice, dtype=float)
    inverse = np.linalg.inv(lattice)
    source = np.asarray(compiled.faces, dtype=float)
    ids = np.asarray(compiled.material_ids, dtype=int)
    fractional = source.copy()
    fractional[:, :, :2] = source[:, :, :2] @ inverse
    keep = np.ones(len(source), dtype=bool)
    for axis in (0, 1):
        for position in (-.5, .5):
            keep &= ~np.all(abs(fractional[:, :, axis] - position) < 1e-10, axis=1)
    limits = ((-nx / 2, nx / 2), (-ny / 2, ny / 2))
    polygons, materials = [], []
    ranges = [range(math.ceil(-count / 2 - .5 + 1e-10), math.floor(count / 2 + .5 - 1e-10) + 1) for count in (nx, ny)]
    for u in ranges[0]:
        for v in ranges[1]:
            for face, material in zip(fractional[keep], ids[keep]):
                polygon = face + [u, v, 0]
                for axis in (0, 1):
                    polygon = _clip(polygon, axis, limits[axis][0], True)
                    polygon = _clip(polygon, axis, limits[axis][1], False)
                if len(polygon) >= 3 and _area(polygon) > 1e-16:
                    polygons.append(polygon)
                    materials.append(int(material))
    bottom = float(np.min(source[:, :, 2]))
    base = float(np.max(source[ids == 1, :, 2]))
    film_level = float(np.max(source[ids == 2, :, 2])) if np.any(ids == 2) else base
    cap_faces, cap_ids = _caps(polygons, materials, limits, base, bottom, film_level)
    polygons.extend(cap_faces)
    materials.extend(cap_ids)
    spacing = max(float(np.min(np.linalg.norm(lattice, axis=1))) / 12, 1e-12)
    result, result_ids = [], []
    for polygon, material in zip(polygons, materials):
        physical = np.asarray(polygon).copy()
        physical[:, :2] = physical[:, :2] @ lattice
        for face in _quads(physical):
            pieces = _subdivide(face, spacing)
            result.append(pieces)
            result_ids.extend([material] * len(pieces))
    faces = np.concatenate(result)
    bounds = np.column_stack((faces.min(axis=(0, 1)), faces.max(axis=(0, 1))))
    return faces, np.asarray(result_ids, dtype=int), compiled.materials, bounds
