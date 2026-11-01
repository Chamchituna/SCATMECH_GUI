from copy import deepcopy
import math

import numpy as np


def _active(document):
    return document.dimension == 1 and document.one_d_mode == "regions"


def _regions(document):
    return {region["id"]: region for region in document.regions if isinstance(region, dict) and isinstance(region.get("id"), str)}


def _fraction(point, a, b, axis=None):
    delta = b - a
    component = {"x": 0, "z": 1}.get(axis)
    if component is not None and abs(delta[component]) > np.linalg.norm(delta) * 1e-12:
        value = (point[component] - a[component]) / delta[component]
    elif component is not None:
        return None
    else:
        length = float(np.dot(delta, delta))
        if length <= 0:
            return None
        value = float(np.dot(point - a, delta)) / length
    return float(np.clip(value, 1e-6, 1 - 1e-6))


def _order(regions):
    ordered, active = [], set()

    def visit(identity):
        if identity in active:
            raise ValueError("Region attachment dependencies must not form a cycle")
        if identity in ordered:
            return
        active.add(identity)
        for attachment in regions[identity].get("attachments", []):
            visit(attachment["target"])
        active.remove(identity)
        ordered.append(identity)

    for identity in regions:
        visit(identity)
    return ordered


def _ancestors(regions, identity):
    result = set()
    pending = [identity]
    while pending:
        for attachment in regions[pending.pop()].get("attachments", []):
            target = attachment["target"]
            if target not in result:
                result.add(target)
                pending.append(target)
    return result


def ensure_spacer_constraints(document):
    if not _active(document) or not isinstance(document.regions, list):
        return document
    if not isinstance(document.period_um, (int, float)) or isinstance(document.period_um, bool) or not math.isfinite(document.period_um) or document.period_um <= 0:
        return document
    regions = _regions(document)
    if not all(name in regions for name in ("film", "core", "left-spacer", "right-spacer")):
        return document
    core, left, right = (regions[name] for name in ("core", "left-spacer", "right-spacer"))
    film = regions["film"]
    try:
        points = [np.asarray(region.get("vertices", []), dtype=float) for region in (core, left, right)]
    except (ValueError, TypeError):
        return document
    if [point.shape for point in points] != [(8, 2), (5, 2), (5, 2)] or not all(np.isfinite(point).all() for point in points):
        return document
    try:
        film_points = np.asarray(film.get("vertices", []), dtype=float)
    except (ValueError, TypeError):
        return document
    if film_points.shape != (4, 2) or not np.isfinite(film_points).all():
        return document
    changed = False
    for spacer, first, second, contact, edge in ((left, 1, 2, 3, [7, 6]), (right, 3, 2, 1, [2, 3])):
        if "attachments" in spacer:
            continue
        a, b = np.asarray(core["vertices"])[edge]
        fraction = _fraction(np.asarray(spacer["vertices"][contact]), a, b)
        spacer["attachments"] = [
            {"vertex": first, "target": "core", "target_vertex": 0 if spacer is left else 1},
            {"vertex": second, "target": "core", "target_vertex": edge[0]},
            {"vertex": contact, "target": "core", "edge": edge, "fraction": fraction if fraction is not None else 4 / 7},
        ]
        changed = True
    if "base_layer" not in film:
        film["base_layer"] = {"thickness": max(float(film_points[:, 1].max()), float(document.period_um) * 1e-6)}
        changed = True
    for region, vertices in ((core, (0, 1)), (left, (0,)), (right, (4,))):
        if not isinstance(region.get("attachments", []), list):
            continue
        attachments = region.setdefault("attachments", [])
        existing = {item.get("vertex") for item in attachments if isinstance(item, dict)}
        for vertex in vertices:
            if vertex not in existing:
                fraction = float(np.clip(region["vertices"][vertex][0] / document.period_um + .5, 1e-6, 1 - 1e-6))
                attachments.append({"vertex": vertex, "target": "film", "edge": [3, 2], "fraction": fraction})
                changed = True
    if changed:
        enforce_region_constraints(document)
    return document


def _index(value, count):
    return isinstance(value, int) and not isinstance(value, bool) and 0 <= value < count


def constraint_errors(document, check_positions=True):
    if not _active(document):
        return []
    errors = []
    regions = _regions(document)
    tolerance = max(abs(float(document.period_um)), 1e-6) * 1e-9
    for region in regions.values():
        attachments = region.get("attachments", [])
        label = str(region.get("name", region["id"]))
        base = region.get("base_layer")
        if base is not None:
            height = base.get("thickness") if isinstance(base, dict) else None
            if not isinstance(height, (int, float)) or isinstance(height, bool) or not math.isfinite(height) or height <= 0 or len(region.get("vertices", [])) != 4:
                errors.append(label + " must be a rectangular substrate film with positive thickness")
            elif check_positions:
                period = document.period_um
                expected = [[-period / 2, 0], [period / 2, 0], [period / 2, height], [-period / 2, height]]
                if not np.allclose(region["vertices"], expected, rtol=0, atol=tolerance):
                    errors.append(label + " must span one period and remain attached to the substrate")
            if attachments:
                errors.append(label + " substrate film cannot depend on another region")
        if not isinstance(attachments, list):
            errors.append(label + " attachments must be a list")
            continue
        seen = set()
        for attachment in attachments:
            prefix = label + " attachment: "
            if not isinstance(attachment, dict) or not _index(attachment.get("vertex"), len(region.get("vertices", []))):
                errors.append(prefix + "invalid source vertex")
                continue
            vertex = attachment["vertex"]
            if vertex in seen:
                errors.append(prefix + "a vertex cannot have multiple attachments")
            seen.add(vertex)
            target_id = attachment.get("target")
            target = regions.get(target_id) if isinstance(target_id, str) else None
            if target is None or target is region:
                errors.append(prefix + "target must be another existing region")
                continue
            points = target.get("vertices", [])
            if "target_vertex" in attachment and "edge" not in attachment:
                target_vertex = attachment["target_vertex"]
                if not _index(target_vertex, len(points)):
                    errors.append(prefix + "invalid target vertex")
                    continue
                expected = points[target_vertex]
            elif "edge" in attachment and "target_vertex" not in attachment:
                edge = attachment["edge"]
                fraction = attachment.get("fraction")
                if not isinstance(edge, list) or len(edge) != 2 or not all(_index(index, len(points)) for index in edge):
                    errors.append(prefix + "invalid target edge")
                    continue
                if (edge[0] - edge[1]) % len(points) not in (1, len(points) - 1):
                    errors.append(prefix + "target edge endpoints must be adjacent")
                    continue
                if not isinstance(fraction, (int, float)) or isinstance(fraction, bool) or not math.isfinite(fraction) or not 0 < fraction < 1:
                    errors.append(prefix + "edge position must be strictly between its endpoints")
                    continue
                a, b = np.asarray(points, dtype=float)[edge]
                if check_positions and np.linalg.norm(b - a) <= tolerance:
                    errors.append(prefix + "target edge must have nonzero length")
                    continue
                expected = a + fraction * (b - a)
            else:
                errors.append(prefix + "specify a target vertex or edge")
                continue
            if check_positions and not np.allclose(region["vertices"][vertex], expected, rtol=0, atol=tolerance):
                errors.append(prefix + "vertex must remain connected to " + str(target.get("name", target_id)))
    if not errors:
        try:
            _order(regions)
        except ValueError as error:
            errors.append(str(error))
    return errors


def enforce_region_constraints(document, previous=None, edited=None, axis=None):
    if not _active(document):
        return document
    errors = constraint_errors(document, check_positions=False)
    if errors:
        raise ValueError("; ".join(errors))
    regions = _regions(document)
    old = _regions(previous) if previous is not None and _active(previous) else {}
    ordered = _order(regions)
    affected = set()
    contact = edited
    while contact is not None and contact[0] in regions:
        affected.add(contact)
        source = regions[contact[0]]
        attachment = next((item for item in source.get("attachments", []) if item["vertex"] == contact[1] and "target_vertex" in item), None)
        if attachment is None:
            break
        point = deepcopy(source["vertices"][contact[1]])
        contact = (attachment["target"], attachment["target_vertex"])
        regions[contact[0]]["vertices"][contact[1]] = point
    for identity in ordered:
        region = regions[identity]
        if "base_layer" not in region:
            continue
        height = region["base_layer"]["thickness"]
        before = old.get(identity)
        if edited is not None and edited[0] == identity:
            if edited[1] in (2, 3):
                height = max(float(region["vertices"][edited[1]][1]), document.period_um * 1e-6)
        elif before is not None and before.get("base_layer") == region["base_layer"] and before["vertices"] != region["vertices"]:
            proposed = max(float(np.ptp(np.asarray(region["vertices"])[:, 1])), document.period_um * 1e-6)
            if not math.isclose(proposed, height, rel_tol=1e-12, abs_tol=document.period_um * 1e-12):
                height = proposed
        region["base_layer"]["thickness"] = height
        period = document.period_um
        region["vertices"] = [[-period / 2, 0.], [period / 2, 0.], [period / 2, height], [-period / 2, height]]
        if before is not None and "base_layer" in before:
            delta = height - before["base_layer"]["thickness"]
            if delta:
                for child in ordered:
                    if identity in _ancestors(regions, child) and child in old and regions[child]["vertices"] == old[child]["vertices"]:
                        regions[child]["vertices"] = (np.asarray(regions[child]["vertices"]) + [0., delta]).tolist()
    for identity in ordered:
        region = regions[identity]
        for attachment in region.get("attachments", []):
            vertex = attachment["vertex"]
            target = regions[attachment["target"]]
            if "target_vertex" in attachment:
                point = target["vertices"][attachment["target_vertex"]]
            else:
                a, b = np.asarray(target["vertices"], dtype=float)[attachment["edge"]]
                before = old.get(identity)
                before_target = old.get(attachment["target"])
                source_moved = before is not None and len(before["vertices"]) == len(region["vertices"]) and not np.array_equal(before["vertices"][vertex], region["vertices"][vertex])
                target_moved = before_target is not None and before_target["vertices"] != target["vertices"]
                if (identity, vertex) in affected or source_moved and not target_moved:
                    fraction = _fraction(np.asarray(region["vertices"][vertex]), a, b, axis if (identity, vertex) in affected else None)
                    if fraction is not None:
                        attachment["fraction"] = fraction
                point = a + attachment["fraction"] * (b - a)
            region["vertices"][vertex] = [float(value) for value in point]
    return document


def _protected(document):
    vertices, edges = {}, {}
    regions = _regions(document)
    for identity, region in regions.items():
        if "base_layer" in region:
            vertices.setdefault(identity, set()).update(range(4))
            edges.setdefault(identity, set()).update(frozenset((index, (index + 1) % 4)) for index in range(4))
        own = set()
        for attachment in region.get("attachments", []):
            own.add(attachment["vertex"])
            target = attachment["target"]
            if "target_vertex" in attachment:
                vertices.setdefault(target, set()).add(attachment["target_vertex"])
            else:
                pair = frozenset(attachment["edge"])
                edges.setdefault(target, set()).add(pair)
                vertices.setdefault(target, set()).update(pair)
        vertices.setdefault(identity, set()).update(own)
        for index in own:
            following = (index + 1) % len(region["vertices"])
            if following in own:
                edges.setdefault(identity, set()).add(frozenset((index, following)))
        vertex_attachments = [item for item in region.get("attachments", []) if "target_vertex" in item]
        for first in vertex_attachments:
            for second in vertex_attachments:
                if first["target"] == second["target"] and (first["vertex"] + 1) % len(region["vertices"]) == second["vertex"]:
                    edges.setdefault(first["target"], set()).add(frozenset((first["target_vertex"], second["target_vertex"])))
    return vertices, edges


def region_vertex_edit_error(document, region_index, vertex_index, operation):
    if not _active(document):
        return ""
    region = document.regions[region_index]
    vertices, edges = _protected(document)
    if operation == "delete" and vertex_index in vertices.get(region["id"], set()):
        return "This vertex maintains a required material connection and cannot be deleted."
    following = (vertex_index + 1) % len(region["vertices"])
    if operation in ("insert", "add") and frozenset((vertex_index, following)) in edges.get(region["id"], set()):
        return "This edge maintains a required material connection. Add vertices on a free edge instead."
    return ""


def region_delete_error(document, region_index):
    region = document.regions[region_index]
    if any(attachment["target"] == region["id"] for item in document.regions for attachment in item.get("attachments", [])):
        return "This region anchors other regions. Delete its attached regions first."
    return ""


def remap_region_vertices(document, region_id, index, delta):
    def remap(value):
        return value + delta if value >= index else value

    for region in document.regions:
        for attachment in region.get("attachments", []):
            if region["id"] == region_id:
                attachment["vertex"] = remap(attachment["vertex"])
            if attachment["target"] == region_id:
                if "target_vertex" in attachment:
                    attachment["target_vertex"] = remap(attachment["target_vertex"])
                else:
                    attachment["edge"] = [remap(value) for value in attachment["edge"]]
