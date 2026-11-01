from __future__ import annotations

import json
import math
from dataclasses import replace
from pathlib import Path

import numpy as np

from geometry_types import GeometryLayer, GeometryRequest, GeometryScene
from scatmech_gratings import CROSS_GRATING_SPECS, ONE_D_GRATING_SPECS, build_default_cross_grating, validate_cross_grating, validate_one_d_grating


def _complex(values):
    array = np.asarray(values)
    if np.iscomplexobj(array):
        return np.asarray(array, dtype=np.complex128)
    if array.shape[-1:] != (2,):
        raise ValueError("Geometry tensors must contain real/imaginary pairs")
    return np.asarray(array[..., 0], dtype=float) + 1j * np.asarray(array[..., 1], dtype=float)


def _lattice(d1, d2, zeta):
    angle = math.radians(zeta)
    result = np.array([[d1, 0.0], [d2 * math.sin(angle), d2 * math.cos(angle)]])
    if not np.isfinite(result).all() or abs(np.linalg.det(result)) < 1e-14:
        raise ValueError("The evaluated lattice is singular")
    return result


def _radial_geometry(raw):
    if raw is None:
        return None
    angles = np.asarray(raw["angles"], dtype=float)
    radii = np.asarray(raw["radii"], dtype=float)
    inside = _complex(raw["inside"])
    outside = _complex(raw["outside"])
    if angles.ndim != 1 or len(angles) < 3 or radii.shape != angles.shape:
        raise ValueError("Native exporter returned inconsistent radial geometry")
    if (not np.isfinite(angles).all() or not np.isfinite(radii).all()
            or np.any(np.diff(angles) <= 0) or angles[-1] - angles[0] >= 2 * math.pi + 1e-12):
        raise ValueError("Native exporter returned invalid radial boundaries")
    if inside.shape != (6,) or outside.shape != (6,) or not np.isfinite(inside).all() or not np.isfinite(outside).all():
        raise ValueError("Native exporter returned invalid radial material tensors")
    return {"angles": angles, "radii": radii, "inside": inside, "outside": outside}


def _leaf_scene(data, wavelength):
    lattice = np.asarray(data["lattice"], dtype=float)
    if lattice.shape != (2, 2) or not np.isfinite(lattice).all() or abs(np.linalg.det(lattice)) < 1e-14:
        raise ValueError("Native exporter returned an invalid lattice")
    layers = []
    z_top = 0.0
    for index, raw in enumerate(data["layers"]):
        thickness = float(raw["thickness"])
        if not np.isfinite(thickness) or thickness < 0:
            raise ValueError("Native exporter returned an invalid layer thickness")
        u_edges = np.asarray(raw["u_edges"], dtype=float)
        v_edges = np.asarray(raw["v_edges"], dtype=float)
        tensors = _complex(raw["tensors"])
        if tensors.shape != (len(v_edges) - 1, len(u_edges) - 1, 6):
            raise ValueError("Native exporter returned inconsistent layer dimensions")
        if not np.isfinite(tensors).all() or any(not np.isfinite(edges).all() or np.any(np.diff(edges) <= 0) for edges in (u_edges, v_edges)):
            raise ValueError("Native exporter returned invalid material samples")
        if thickness:
            layers.append(GeometryLayer(z_top, z_top - thickness, u_edges, v_edges, tensors, lattice.copy(), label=f'{data["model"]} layer {index + 1}', continuous=bool(raw.get("continuous", False)), radial=_radial_geometry(raw.get("radial"))))
        z_top -= thickness
    medium_i = _complex(data["medium_i"])
    medium_t = _complex(data["medium_t"])
    if medium_i.shape != (6,) or medium_t.shape != (6,) or not np.isfinite(medium_i).all() or not np.isfinite(medium_t).all():
        raise ValueError("Native exporter returned invalid incident or substrate material")
    return GeometryScene(data["model"], wavelength, lattice, layers, medium_i, medium_t, int(data["dimension"]), files=list(data.get("files", [])), notes=list(data.get("notes", [])))


def _transform_one_d(scene, normal, tangent):
    transform = np.array([normal, tangent], dtype=float)
    return [replace(layer, basis=layer.basis @ transform, offset=layer.offset @ transform) for layer in scene.layers]


def _same_tensor(first, second, message):
    if not np.allclose(first, second, rtol=1e-9, atol=1e-12):
        raise ValueError(message)


def _tensor_notes(layers, notes):
    if any(np.any(layer.tensors[..., :3] != layer.tensors[..., :1]) or np.any(layer.tensors[..., 3:] != layer.tensors[..., 3:4]) for layer in layers):
        return list(dict.fromkeys(notes + ["Anisotropic ε/μ components retain their native material axes; geometry rotations are recorded in each layer basis."]))
    return notes


def build_scene(request: GeometryRequest, progress=None) -> GeometryScene:
    from geometry_native import export_leaf

    node = request.grating
    model = node.get("model") if isinstance(node, dict) else None
    validate = validate_one_d_grating if model in ONE_D_GRATING_SPECS else validate_cross_grating
    errors = validate(node)
    if errors:
        raise ValueError("\n".join(errors))
    if not np.isfinite(request.wavelength_um) or request.wavelength_um <= 0:
        raise ValueError("Preview wavelength must be positive and finite")
    if not np.isfinite(request.incident_deg) or not np.isfinite(request.rotation_deg):
        raise ValueError("Incident and rotation angles must be finite")

    def evaluate(current):
        model = current["model"]
        spec = (ONE_D_GRATING_SPECS if model in ONE_D_GRATING_SPECS else CROSS_GRATING_SPECS)[model]
        params = {field["name"]: field["default"] for field in spec["fields"] if field["kind"] == "param"}
        params.update(current.get("params", {}))
        children = current.get("children", {})
        if progress:
            progress(f"Evaluating {model}")
        if model in ONE_D_GRATING_SPECS or model not in {"OneD_CrossGrating", "Overlaid_CrossGrating", "Overlaid_1D_CrossGrating"}:
            return _leaf_scene(export_leaf(current, request.wavelength_um, request.resolution, progress), request.wavelength_um)
        interface = build_default_cross_grating("Null_CrossGrating")
        interface["params"].update(medium_i=params["medium_i"], medium_t=params["medium_t"])
        media = _leaf_scene(export_leaf(interface, request.wavelength_um, request.resolution, progress), request.wavelength_um)
        if model == "OneD_CrossGrating":
            child = evaluate(children["grating"])
            _same_tensor(media.medium_i, child.medium_i, "Evaluated incident material must match the child grating")
            _same_tensor(media.medium_t, child.medium_t, "Evaluated transmission material must match the child grating")
            zeta = math.radians(float(params["zeta"]))
            cosine, sine = math.cos(zeta), math.sin(zeta)
            if abs(cosine) < 1e-12:
                raise ValueError("OneD_CrossGrating requires a nonsingular zeta")
            lattice = _lattice(child.lattice[0, 0] / cosine, float(params["d2"]), float(params["zeta"]))
            layers = _transform_one_d(child, [cosine, -sine], [sine, cosine])
            return GeometryScene(model, request.wavelength_um, lattice, layers, child.medium_i, child.medium_t, 2, files=list(dict.fromkeys(media.files + child.files)), notes=child.notes)
        top = evaluate(children["top"])
        bottom = evaluate(children["bottom"])
        _same_tensor(media.medium_i, top.medium_i, "Evaluated incident material must match the top grating")
        _same_tensor(media.medium_t, bottom.medium_t, "Evaluated transmission material must match the bottom grating")
        _same_tensor(top.medium_t, bottom.medium_i, "Evaluated top transmission material must match bottom incident material")
        separation = float(params["separation"])
        if separation < 0:
            raise ValueError("Overlay separation must be nonnegative")
        if model == "Overlaid_CrossGrating":
            if not np.allclose(top.lattice, bottom.lattice, rtol=1e-9, atol=1e-12):
                raise ValueError("Evaluated top and bottom lattices differ; custom-file lattice values override form parameters")
            lattice = bottom.lattice.copy()
            d2 = np.linalg.norm(lattice[1])
            offset = np.array([float(params["overlay1"]), 0.0]) + float(params["overlay2"]) * lattice[1] / d2
            top_layers = [replace(layer, offset=layer.offset + offset) for layer in top.layers]
            bottom_layers = bottom.layers
        else:
            angle = math.radians(float(params["angle"]))
            sine, cosine = math.sin(angle), math.cos(angle)
            if abs(sine) < 1e-12:
                raise ValueError("Overlay angle produces a singular lattice")
            lattice = _lattice(top.lattice[0, 0] / sine, bottom.lattice[0, 0] / sine, 90.0 - float(params["angle"]))
            top_layers = _transform_one_d(top, [sine, -cosine], [cosine, sine])
            bottom_layers = _transform_one_d(bottom, [0.0, 1.0], [-1.0, 0.0])
        layers = list(top_layers)
        if separation:
            layers.append(GeometryLayer(-top.thickness, -top.thickness - separation, np.array([-.5, .5]), np.array([-.5, .5]), top.medium_t.reshape(1, 1, 6), lattice.copy(), label="Overlay separation"))
        shift = top.thickness + separation
        layers.extend(replace(layer, z_top=layer.z_top - shift, z_bottom=layer.z_bottom - shift) for layer in bottom_layers)
        return GeometryScene(model, request.wavelength_um, lattice, layers, top.medium_i, bottom.medium_t, 2, files=list(dict.fromkeys(media.files + top.files + bottom.files)), notes=list(dict.fromkeys(top.notes + bottom.notes)))

    scene = evaluate(node)
    scene.incident_deg = float(request.incident_deg)
    scene.rotation_deg = float(request.rotation_deg)
    scene.configuration = int(request.configuration)
    scene.notes = _tensor_notes(scene.layers, scene.notes)
    return scene


def save_scene(scene: GeometryScene, path):
    metadata = {"version": 1, "model": scene.model, "wavelength_um": scene.wavelength_um, "dimension": scene.dimension, "incident_deg": scene.incident_deg, "rotation_deg": scene.rotation_deg, "configuration": scene.configuration, "files": scene.files, "notes": scene.notes, "layers": []}
    arrays = {"lattice": scene.lattice, "medium_i": scene.medium_i, "medium_t": scene.medium_t}
    for index, layer in enumerate(scene.layers):
        metadata["layers"].append({"z_top": layer.z_top, "z_bottom": layer.z_bottom, "label": layer.label, "continuous": layer.continuous, "has_radial": layer.radial is not None})
        for name in ("u_edges", "v_edges", "tensors", "basis", "offset"):
            arrays[f"layer_{index}_{name}"] = getattr(layer, name)
        if layer.radial is not None:
            for name, value in layer.radial.items():
                arrays[f"layer_{index}_radial_{name}"] = value
    arrays["metadata"] = np.asarray(json.dumps(metadata))
    with Path(path).open("wb") as handle:
        np.savez_compressed(handle, **arrays)


def load_scene(path) -> GeometryScene:
    with np.load(path, allow_pickle=False) as archive:
        metadata = json.loads(str(archive["metadata"].item()))
        if metadata.pop("version") != 1:
            raise ValueError("Unsupported geometry scene version")
        layers = []
        for index, raw in enumerate(metadata.pop("layers")):
            radial = None
            if raw.pop("has_radial", False):
                radial = _radial_geometry({name: archive[f"layer_{index}_radial_{name}"].copy()
                                           for name in ("angles", "radii", "inside", "outside")})
            layers.append(GeometryLayer(**raw, radial=radial, **{
                name: archive[f"layer_{index}_{name}"].copy()
                for name in ("u_edges", "v_edges", "tensors", "basis", "offset")
            }))
        return GeometryScene(**metadata, layers=layers, lattice=archive["lattice"].copy(), medium_i=archive["medium_i"].copy(), medium_t=archive["medium_t"].copy())
