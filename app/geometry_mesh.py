from __future__ import annotations

import math

import numpy as np

from geometry_types import GeometryMesh


def _radial_values(radial, angles):
    source = np.asarray(radial["angles"])
    if source[-1] - source[0] >= 2 * np.pi - 1e-10:
        return np.interp(angles, source, radial["radii"])
    return np.interp(angles, source, radial["radii"], period=2 * np.pi)


def sample_layer(layer, points, smooth_boundaries=False):
    local = (np.asarray(points) - layer.offset) @ np.linalg.inv(layer.basis)
    u0, u1 = layer.u_edges[0], layer.u_edges[-1]
    v0, v1 = layer.v_edges[0], layer.v_edges[-1]
    u = (local[..., 0] - u0) % (u1 - u0) + u0
    v = (local[..., 1] - v0) % (v1 - v0) + v0
    if smooth_boundaries and getattr(layer, "radial", None) is not None:
        xy = np.stack([u, v], axis=-1) @ layer.basis
        angles = np.arctan2(xy[..., 1], xy[..., 0])
        radius = _radial_values(layer.radial, angles)
        return np.where((np.linalg.norm(xy, axis=-1) <= radius)[..., None], layer.radial["inside"], layer.radial["outside"])
    ix = np.clip(np.searchsorted(layer.u_edges, u, side="right") - 1, 0, len(layer.u_edges) - 2)
    iy = np.clip(np.searchsorted(layer.v_edges, v, side="right") - 1, 0, len(layer.v_edges) - 2)
    return layer.tensors[iy, ix]


def _grid_edges(scene, layer, cap):
    mapping = scene.lattice @ np.linalg.inv(layer.basis)
    offset = layer.offset @ np.linalg.inv(layer.basis)
    result = []
    simplified = False
    for axis in (0, 1):
        boundaries = [-.5, .5]
        aligned = False
        for source_axis, edges in enumerate((layer.u_edges, layer.v_edges)):
            if len(edges) == 2:
                continue
            scale = mapping[axis, source_axis]
            if abs(mapping[1 - axis, source_axis]) < 1e-9 and abs(scale) > 1e-10:
                aligned = True
                period = edges[-1] - edges[0]
                low = min(-.5 * scale - offset[source_axis], .5 * scale - offset[source_axis])
                high = max(-.5 * scale - offset[source_axis], .5 * scale - offset[source_axis])
                for shift in range(math.floor((low - edges[-1]) / period), math.ceil((high - edges[0]) / period) + 1):
                    values = (edges + shift * period + offset[source_axis]) / scale
                    boundaries.extend(values[(values > -.5 + 1e-12) & (values < .5 - 1e-12)].tolist())
        if not aligned:
            needs_sample = any(len(edges) > 2 and abs(mapping[axis, j]) > 1e-9 for j, edges in enumerate((layer.u_edges, layer.v_edges)))
            if needs_sample:
                boundaries = np.linspace(-.5, .5, cap + 1).tolist()
                simplified = True
        values = np.unique(np.round(boundaries, 13))
        if len(values) > cap + 1:
            indices = np.unique(np.rint(np.linspace(0, len(values) - 1, cap + 1)).astype(int))
            values = values[indices]
            simplified = True
        result.append(values)
    return result[0], result[1], simplified


def _polygon_clip(vertices, axis, cutoff):
    if axis is None:
        return list(vertices)
    result = []
    previous = vertices[-1]
    previous_inside = previous[axis] <= cutoff + 1e-12
    for current in vertices:
        inside = current[axis] <= cutoff + 1e-12
        if inside != previous_inside:
            fraction = (cutoff - previous[axis]) / (current[axis] - previous[axis])
            result.append(previous + fraction * (current - previous))
        if inside:
            result.append(current)
        previous, previous_inside = current, inside
    return result


def _quads(polygon):
    if len(polygon) < 3:
        return []
    if len(polygon) == 3:
        return [[polygon[0], polygon[1], polygon[2], polygon[2]]]
    if len(polygon) == 4:
        return [polygon]
    return [[polygon[0], polygon[index], polygon[index + 1], polygon[index + 1]] for index in range(1, len(polygon) - 1)]


def _cap(vertices, axis, cutoff):
    if axis is None or min(vertices[:, axis]) >= cutoff - 1e-12 or max(vertices[:, axis]) <= cutoff + 1e-12:
        return []
    points = []
    for first, second in ((0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4), (0, 4), (1, 5), (2, 6), (3, 7)):
        a, b = vertices[first], vertices[second]
        if (a[axis] - cutoff) * (b[axis] - cutoff) < 0:
            points.append(a + (cutoff - a[axis]) / (b[axis] - a[axis]) * (b - a))
        elif abs(a[axis] - cutoff) < 1e-12:
            points.append(a)
        elif abs(b[axis] - cutoff) < 1e-12:
            points.append(b)
    if len(points) < 3:
        return []
    points = np.unique(np.round(points, 13), axis=0)
    other = [value for value in range(3) if value != axis]
    planar = points[:, other] - points[:, other].mean(axis=0)
    return list(points[np.argsort(np.arctan2(planar[:, 1], planar[:, 0]))])


def _repeat_edges(edges, count):
    low, high = -count / 2, count / 2
    values = np.concatenate([edges + shift for shift in range(math.floor(low - .5), math.ceil(high + .5) + 1)])
    return np.unique(np.round(np.concatenate([[low], values[(values > low + 1e-12) & (values < high - 1e-12)], [high]]), 13))


def _clip_xy(polygon, normal, cutoff):
    if len(polygon) < 3:
        return np.empty((0, 2))
    result = []
    previous = polygon[-1]
    previous_distance = float(previous @ normal - cutoff)
    for current in polygon:
        distance = float(current @ normal - cutoff)
        if (distance <= 1e-12) != (previous_distance <= 1e-12):
            fraction = previous_distance / (previous_distance - distance)
            result.append(previous + fraction * (current - previous))
        if distance <= 1e-12:
            result.append(current)
        previous, previous_distance = current, distance
    if not result:
        return np.empty((0, 2))
    points = np.asarray(result)
    keep = np.linalg.norm(points - np.roll(points, 1, axis=0), axis=1) > 1e-12
    return points[keep]


def _radial_indices(layer, limit):
    angles = np.asarray(layer.radial["angles"])
    if len(angles) <= limit:
        return np.arange(len(angles))
    corners = np.array([[-.5, -.5], [.5, -.5], [.5, .5], [-.5, .5]]) @ layer.basis
    corner_angles = np.arctan2(corners[:, 1], corners[:, 0])
    retained = [int(np.argmin(np.abs(np.angle(np.exp(1j * (angles - angle)))))) for angle in corner_angles]
    retained.extend([0, len(angles) - 1])
    retained.extend(np.rint(np.linspace(0, len(angles), max(8, limit - 6), endpoint=False)).astype(int).tolist())
    return np.unique(retained)


def _aligned_radial(first, second):
    if second is None or getattr(second, "radial", None) is None or not np.allclose(first.basis, second.basis, atol=1e-12, rtol=1e-10):
        return False
    offset = (first.offset - second.offset) @ np.linalg.inv(first.basis)
    return np.allclose(offset, np.rint(offset), atol=1e-10, rtol=0)


def _uniform_layer(layer):
    if layer is None:
        return True
    radial = getattr(layer, "radial", None)
    if radial is not None:
        return np.array_equal(radial["inside"], radial["outside"])
    return bool(np.all(layer.tensors == layer.tensors[0, 0]))


def _radial_segments(layer, limit):
    angles = np.asarray(layer.radial["angles"])[_radial_indices(layer, limit)]
    return len(angles) - int(angles[-1] - angles[0] >= 2 * np.pi - 1e-10)


def _radial_tiles(layer, scene, repeats):
    nx, ny = repeats
    corners = np.array([[-nx / 2, -ny / 2], [nx / 2, -ny / 2], [nx / 2, ny / 2], [-nx / 2, ny / 2]]) @ scene.lattice
    local = (corners - layer.offset) @ np.linalg.inv(layer.basis)
    low = np.ceil(local.min(axis=0) - .5 + 1e-12).astype(int)
    high = np.floor(local.max(axis=0) + .5 - 1e-12).astype(int)
    return (np.array([iu, iv]) @ layer.basis + layer.offset for iv in range(low[1], high[1] + 1) for iu in range(low[0], high[0] + 1))


def _radial_polygons(layer, scene, repeats, limit, adjacent):
    indices = _radial_indices(layer, limit)
    angles = np.asarray(layer.radial["angles"])[indices]
    directions = np.column_stack([np.cos(angles), np.sin(angles)])
    local_directions = directions @ np.linalg.inv(layer.basis)
    outer_radii = .5 / np.max(np.abs(local_directions), axis=1)
    boundaries = [np.zeros(len(angles)), outer_radii]
    for neighbor in [layer] + [item for item in adjacent if _aligned_radial(layer, item)]:
        radii = _radial_values(neighbor.radial, angles)
        boundaries.append(np.clip(radii, 0, outer_radii))
    contour = boundaries[2]
    boundaries = np.sort(np.asarray(boundaries), axis=0)
    for center in _radial_tiles(layer, scene, repeats):
        closed = angles[-1] - angles[0] >= 2 * np.pi - 1e-10
        for index in range(len(angles) - int(closed)):
            following = (index + 1) % len(angles)
            direction = directions[[index, following]]
            for band in range(len(boundaries) - 1):
                inner = boundaries[band, [index, following]]
                outer = boundaries[band + 1, [index, following]]
                if np.max(outer - inner) <= 1e-12:
                    continue
                points = np.array([direction[0] * inner[0], direction[0] * outer[0], direction[1] * outer[1], direction[1] * inner[1]]) + center
                boundary = contour[[index, following]]
                tensor = layer.radial["inside"] if np.mean(inner + outer) / 2 <= np.mean(boundary) else layer.radial["outside"]
                sides = []
                if np.any(layer.radial["inside"] != layer.radial["outside"]):
                    if np.allclose(outer, boundary, rtol=0, atol=1e-12) and np.any(boundary < outer_radii[[index, following]] - 1e-12):
                        sides.append((points[1], points[2]))
                    if np.allclose(inner, boundary, rtol=0, atol=1e-12):
                        sides.append((points[3], points[0]))
                if np.allclose(outer, outer_radii[[index, following]], rtol=0, atol=1e-12):
                    midpoint = (points[1] + points[2]) / 2
                    edge = points[2] - points[1]
                    normal = np.array([-edge[1], edge[0]])
                    normal *= max(np.linalg.norm(layer.basis, axis=1)) * 1e-8 / max(np.linalg.norm(normal), 1e-12)
                    across = sample_layer(layer, np.array([midpoint - normal, midpoint + normal]), True)
                    if np.any(across[0] != across[1]):
                        sides.append((points[1], points[2]))
                yield points, tensor, sides


def _radial_seams(layer, scene, repeats):
    radial = layer.radial
    angles = np.asarray(radial["angles"])
    if angles[-1] - angles[0] < 2 * np.pi - 1e-10 or np.array_equal(radial["inside"], radial["outside"]):
        return
    direction = np.array([np.cos(angles[0]), np.sin(angles[0])])
    extent = .5 / np.max(np.abs(direction @ np.linalg.inv(layer.basis)))
    radii = np.sort(np.clip(np.asarray(radial["radii"])[[0, -1]], 0, extent))
    if radii[1] - radii[0] < 1e-12:
        return
    for center in _radial_tiles(layer, scene, repeats):
        yield center + radii[:, None] * direction


def _clip_line(line, planes):
    first, second = np.asarray(line)
    for normal, cutoff in planes:
        distance = np.array([first @ normal, second @ normal]) - cutoff
        if np.all(distance > 1e-12):
            return None
        if (distance[0] > 1e-12) != (distance[1] > 1e-12):
            intersection = first + distance[0] / (distance[0] - distance[1]) * (second - first)
            if distance[0] > 1e-12:
                first = intersection
            else:
                second = intersection
    return np.array([first, second])


def _on_segment(first, second, start, end):
    edge = end - start
    length = np.linalg.norm(edge)
    if length < 1e-12:
        return False
    difference = np.array([first, second]) - start
    cross = difference[:, 0] * edge[1] - difference[:, 1] * edge[0]
    projection = difference @ edge / (length * length)
    return bool(np.all(np.abs(cross) <= 1e-10 * length) and np.all(projection >= -1e-10) and np.all(projection <= 1 + 1e-10))


def build_mesh(scene, repeats=(3, 3), slice_axis=None, slice_fraction=1.0, max_faces=30000, smooth_boundaries=True) -> GeometryMesh:
    repeats = tuple(int(value) for value in repeats)
    if len(repeats) != 2 or any(value < 1 or value > 25 for value in repeats):
        raise ValueError("Repetitions must be between 1 and 25")
    if max_faces < 12:
        raise ValueError("Mesh face limit must be at least 12")
    if slice_axis not in (None, "x", "y", "z"):
        raise ValueError("Cross-section axis must be x, y, z, or None")
    if not np.isfinite(slice_fraction) or not 0 < slice_fraction <= 1:
        raise ValueError("Cross-section fraction must be greater than zero and at most one")
    nx_repeat, ny_repeat = repeats
    cell_corners = np.array([[-nx_repeat / 2, -ny_repeat / 2], [nx_repeat / 2, -ny_repeat / 2], [nx_repeat / 2, ny_repeat / 2], [-nx_repeat / 2, ny_repeat / 2]]) @ scene.lattice
    slab = max(np.linalg.norm(scene.lattice, axis=1).min() * .15, scene.thickness * .15, .01)
    bounds = np.array([[cell_corners[:, 0].min(), cell_corners[:, 0].max()], [cell_corners[:, 1].min(), cell_corners[:, 1].max()], [-scene.thickness - slab, 0.0]])
    axis = None if slice_axis is None or slice_fraction == 1.0 else "xyz".index(slice_axis)
    cutoff = None if axis is None else bounds[axis, 0] + slice_fraction * np.ptp(bounds[axis])
    materials = [("Incident medium", scene.medium_i.copy()), ("Substrate", scene.medium_t.copy())]
    material_lookup = {np.ascontiguousarray(scene.medium_i).tobytes(): 0}
    if not np.array_equal(scene.medium_i, scene.medium_t):
        material_lookup[np.ascontiguousarray(scene.medium_t).tobytes()] = 1
    continuous_ids = {}
    radial_layers = [layer for layer in scene.layers if smooth_boundaries and getattr(layer, "radial", None) is not None]
    radial_count = len(radial_layers)
    native_angles = max((len(layer.radial["angles"]) for layer in radial_layers), default=0)
    angular_limit = max(12, min(native_angles, int(max_faces / max(6 * radial_count * nx_repeat * ny_repeat, 1))))
    display_planes = []
    inverse_lattice = np.linalg.inv(scene.lattice)
    for coordinate, count in enumerate(repeats):
        display_planes.extend([(inverse_lattice[:, coordinate], count / 2), (-inverse_lattice[:, coordinate], count / 2)])
    if axis in (0, 1):
        display_planes.append((np.eye(2)[axis], cutoff))

    def material_id(tensor):
        key = np.ascontiguousarray(tensor).tobytes()
        if key not in material_lookup:
            material_lookup[key] = len(materials)
            materials.append((f"Material {len(materials) - 1}", tensor.copy()))
        return material_lookup[key]

    cap = min(64, max(1, int(math.sqrt(max_faces / max(2 * len(scene.layers) * nx_repeat * ny_repeat, 1)))))
    source_cells = sum(layer.tensors.shape[0] * layer.tensors.shape[1] for layer in scene.layers)
    attempt = 0
    while True:
        faces, identifiers, face_tensors = [], [], []
        display_cells = 0
        simplified = False

        def emit(vertices, selected_faces, tensor, identifier):
            if axis is not None and np.min(vertices[:, axis]) >= cutoff - 1e-12:
                return
            for indices in selected_faces:
                for quad in _quads(_polygon_clip(vertices[list(indices)], axis, cutoff)):
                    faces.append(quad)
                    identifiers.append(identifier)
                    face_tensors.append(tensor)
            for quad in _quads(_cap(vertices, axis, cutoff)):
                faces.append(quad)
                identifiers.append(identifier)
                face_tensors.append(tensor)

        def prism(xy, bottom, top):
            return np.vstack([np.column_stack([xy, np.full(4, bottom)]), np.column_stack([xy, np.full(4, top)])])

        def emit_radial(layer, above, below):
            if axis == 2 and layer.z_bottom >= cutoff - 1e-12:
                return
            top = min(layer.z_top, cutoff) if axis == 2 else layer.z_top
            bottom = layer.z_bottom
            for polygon, tensor, selected_sides in _radial_polygons(layer, scene, repeats, angular_limit, [above, below]):
                for normal, boundary in display_planes:
                    polygon = _clip_xy(polygon, normal, boundary)
                    if len(polygon) < 3:
                        break
                if len(polygon) < 3:
                    continue
                area = np.sum(polygon[:, 0] * np.roll(polygon[:, 1], -1) - polygon[:, 1] * np.roll(polygon[:, 0], -1))
                if abs(area) < 1e-16:
                    continue
                identifier = material_id(tensor)
                center = polygon.mean(axis=0)
                above_tensor = sample_layer(above, center, smooth_boundaries) if above is not None and abs(above.z_bottom - layer.z_top) < 1e-10 else scene.medium_i
                below_tensor = sample_layer(below, center, smooth_boundaries) if below is not None and abs(below.z_top - layer.z_bottom) < 1e-10 else scene.medium_t

                def append_polygon(points):
                    for quad in _quads(points):
                        faces.append(quad)
                        identifiers.append(identifier)
                        face_tensors.append(tensor)

                if (axis == 2 and abs(top - cutoff) < 1e-12) or np.any(tensor != above_tensor) or not (_uniform_layer(above) or _aligned_radial(layer, above)):
                    append_polygon(np.column_stack([polygon, np.full(len(polygon), top)]))
                if np.any(tensor != below_tensor) or not (_uniform_layer(below) or _aligned_radial(layer, below)):
                    append_polygon(np.column_stack([polygon[::-1], np.full(len(polygon), bottom)]))
                for edge_index, first in enumerate(polygon):
                    second = polygon[(edge_index + 1) % len(polygon)]
                    boundary = any(abs(first @ normal - value) < 1e-10 and abs(second @ normal - value) < 1e-10 for normal, value in display_planes)
                    if boundary or any(_on_segment(first, second, start, end) for start, end in selected_sides):
                        append_polygon(np.array([[first[0], first[1], bottom], [second[0], second[1], bottom], [second[0], second[1], top], [first[0], first[1], top]]))
            for seam in _radial_seams(layer, scene, repeats):
                clipped = _clip_line(seam, display_planes)
                if clipped is None or np.linalg.norm(clipped[1] - clipped[0]) < 1e-12:
                    continue
                first, second = clipped
                quad = np.array([[first[0], first[1], bottom], [second[0], second[1], bottom], [second[0], second[1], top], [first[0], first[1], top]])
                for tensor in (layer.radial["inside"], layer.radial["outside"]):
                    faces.append(quad)
                    identifiers.append(material_id(tensor))
                    face_tensors.append(tensor)

        face_indices = [(0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]
        if scene.layers:
            substrate_u = np.linspace(-nx_repeat / 2, nx_repeat / 2, min(12, cap) * nx_repeat + 1)
            substrate_v = np.linspace(-ny_repeat / 2, ny_repeat / 2, min(12, cap) * ny_repeat + 1)
            for iy in range(len(substrate_v) - 1):
                for ix in range(len(substrate_u) - 1):
                    xy = np.array([[substrate_u[ix], substrate_v[iy]], [substrate_u[ix + 1], substrate_v[iy]], [substrate_u[ix + 1], substrate_v[iy + 1]], [substrate_u[ix], substrate_v[iy + 1]]]) @ scene.lattice
                    selected = face_indices[:2]
                    if iy == 0:
                        selected = selected + [face_indices[2]]
                    if ix == len(substrate_u) - 2:
                        selected = selected + [face_indices[3]]
                    if iy == len(substrate_v) - 2:
                        selected = selected + [face_indices[4]]
                    if ix == 0:
                        selected = selected + [face_indices[5]]
                    emit(prism(xy, bounds[2, 0], -scene.thickness), selected, scene.medium_t, 1)
        else:
            emit(prism(cell_corners, bounds[2, 0], -scene.thickness), face_indices, scene.medium_t, 1)
        grids = [_grid_edges(scene, layer, cap) for layer in scene.layers]
        for layer_index, layer in enumerate(scene.layers):
            above = scene.layers[layer_index - 1] if layer_index else None
            below = scene.layers[layer_index + 1] if layer_index + 1 < len(scene.layers) else None
            if smooth_boundaries and getattr(layer, "radial", None) is not None:
                display_cells += _radial_segments(layer, angular_limit)
                emit_radial(layer, above, below)
                continue
            neighbors = grids[max(0, layer_index - 1):min(len(grids), layer_index + 2)]
            reduced = any(grid[2] for grid in neighbors)
            axes = []
            for grid_axis in (0, 1):
                edges = np.unique(np.concatenate([grid[grid_axis] for grid in neighbors]))
                if len(edges) > cap + 1:
                    edges = edges[np.unique(np.rint(np.linspace(0, len(edges) - 1, cap + 1)).astype(int))]
                    reduced = True
                if len(edges) == 2:
                    edges = np.linspace(edges[0], edges[-1], min(12, cap) + 1)
                axes.append(edges)
            u, v = axes
            simplified |= reduced
            display_cells += (len(u) - 1) * (len(v) - 1)
            if layer.continuous and layer_index not in continuous_ids:
                continuous_ids[layer_index] = len(materials)
                materials.append((f"Continuous index: {layer.label}", layer.tensors.reshape(-1, 6)[len(layer.tensors.reshape(-1, 6)) // 2].copy()))
            us = _repeat_edges(u, nx_repeat)
            vs = _repeat_edges(v, ny_repeat)
            center_u, center_v = np.meshgrid((us[:-1] + us[1:]) / 2, (vs[:-1] + vs[1:]) / 2)
            points = np.stack([center_u, center_v], axis=-1) @ scene.lattice
            tensors = sample_layer(layer, points, smooth_boundaries)
            above_tensors = sample_layer(above, points, smooth_boundaries) if above is not None and abs(above.z_bottom - layer.z_top) < 1e-10 else np.broadcast_to(scene.medium_i, tensors.shape)
            below_tensors = sample_layer(below, points, smooth_boundaries) if below is not None and abs(below.z_top - layer.z_bottom) < 1e-10 else np.broadcast_to(scene.medium_t, tensors.shape)
            different_above = np.any(tensors != above_tensors, axis=-1)
            different_below = np.any(tensors != below_tensors, axis=-1)
            if smooth_boundaries and above is not None and getattr(above, "radial", None) is not None and abs(above.z_bottom - layer.z_top) < 1e-10:
                different_above[:] = True
            if smooth_boundaries and below is not None and getattr(below, "radial", None) is not None and abs(below.z_top - layer.z_bottom) < 1e-10:
                different_below[:] = True
            for iy in range(len(vs) - 1):
                for ix in range(len(us) - 1):
                    tensor = tensors[iy, ix]
                    identifier = continuous_ids[layer_index] if layer.continuous else material_id(tensor)
                    xy = np.array([[us[ix], vs[iy]], [us[ix + 1], vs[iy]], [us[ix + 1], vs[iy + 1]], [us[ix], vs[iy + 1]]]) @ scene.lattice
                    selected = []
                    if different_below[iy, ix] or layer.continuous:
                        selected.append(face_indices[0])
                    if different_above[iy, ix] or layer.continuous or (axis == 2 and abs(layer.z_top - cutoff) < 1e-12):
                        selected.append(face_indices[1])
                    for face_index, jy, jx in ((2, iy - 1, ix), (3, iy, ix + 1), (4, iy + 1, ix), (5, iy, ix - 1)):
                        exterior = jy < 0 or jx < 0 or jy >= tensors.shape[0] or jx >= tensors.shape[1]
                        if exterior or (not layer.continuous and np.any(tensor != tensors[jy, jx])):
                            selected.append(face_indices[face_index])
                    emit(prism(xy, layer.z_bottom, layer.z_top), selected, tensor, identifier)
        if len(faces) <= max_faces:
            break
        attempt += 1
        if cap == 1 and (not radial_count or angular_limit == 12):
            raise ValueError(f"{len(scene.layers)} layers and {nx_repeat} × {ny_repeat} repetitions require more than {max_faces:,} faces. Reduce repetitions or increase the display face limit.")
        cap = max(1, min(cap - 1, int(cap * math.sqrt(max_faces / len(faces)) * .85)))
        if radial_count:
            angular_limit = max(12, min(angular_limit - 1, int(angular_limit * max_faces / len(faces) * .85)))
    lattice_lines = []
    for iu in _repeat_edges(np.array([-.5, .5]), nx_repeat):
        first = np.array([iu, -ny_repeat / 2]) @ scene.lattice
        last = np.array([iu, ny_repeat / 2]) @ scene.lattice
        lattice_lines.append(np.column_stack([np.vstack([first, last]), [0.0, 0.0]]))
    for iv in _repeat_edges(np.array([-.5, .5]), ny_repeat):
        first = np.array([-nx_repeat / 2, iv]) @ scene.lattice
        last = np.array([nx_repeat / 2, iv]) @ scene.lattice
        lattice_lines.append(np.column_stack([np.vstack([first, last]), [0.0, 0.0]]))
    if axis is not None:
        bounds[axis, 1] = cutoff
        lattice_lines = [np.asarray(line) for original in lattice_lines if len(line := _polygon_clip(original, axis, cutoff)) >= 2]
    notes = list(scene.notes)
    notes.append(f"Display: {display_cells:,} cells per repeated stack; exported samples: {source_cells:,} cells; {len(faces):,} faces.")
    if radial_count:
        displayed = [_radial_segments(layer, angular_limit) for layer in radial_layers]
        native = [_radial_segments(layer, len(layer.radial["angles"])) for layer in radial_layers]
        notes.append(f"Smooth cylinder boundaries: {min(displayed)}–{max(displayed)} contour segments per native layer from SCATMECH-evaluated radii. Solver grids and layer steps are unchanged.")
        if any(first < second for first, second in zip(displayed, native)):
            notes.append("Display contour resolution reduced to fit the face limit; fine angular features may be omitted.")
    if simplified or (attempt and radial_count < len(scene.layers)):
        notes.append("Display mesh simplified by sampling exported cells; narrow features may be omitted. Solver grid settings and layer thicknesses are unchanged.")
    notes.append(f"Substrate slab depth {slab:g} µm is a finite illustration of the semi-infinite medium.")
    if scene.dimension == 1:
        notes.append(f"Invariant direction shown over a finite length of {np.ptp(bounds[1]):g} µm.")
    return GeometryMesh(np.asarray(faces, dtype=float).reshape(-1, 4, 3), np.asarray(identifiers, dtype=int), np.asarray(face_tensors, dtype=complex).reshape(-1, 6), materials, bounds, lattice_lines, notes)
