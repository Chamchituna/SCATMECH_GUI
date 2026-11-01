from __future__ import annotations

import copy
import json
import os
from pathlib import Path
import sys

import numpy as np


def main(argv=None):
    args = sys.argv[1:] if argv is None else argv
    if len(args) != 2:
        print("Expected profile request and output directory.", file=sys.stderr)
        return 2
    if os.name == "posix":
        try:
            os.setsid()
        except OSError:
            if os.getpgrp() != os.getpid():
                return 2
    try:
        from custom_profile import compile_profile, native_material_value, profile_from_dict, remap_boundary_media, validate_profile
        from geometry_native import referenced_files
        from geometry_scene import build_scene, save_scene
        from geometry_types import GeometryRequest

        payload = json.loads(Path(args[0]).read_text(encoding="utf-8"))
        folder = Path(args[1])
        os.environ["SCATMECH_GEOMETRY_WORKDIR"] = str(folder)
        document = profile_from_dict(payload["document"])
        errors = validate_profile(document)
        if errors:
            raise ValueError("\n".join(errors))
        compiled = compile_profile(document)
        context = payload["context"]
        containing = None
        input_files = set(compiled.files)
        profiles = [(document, compiled.fingerprint)]
        if context.get("validation_grating") is not None:
            containing = copy.deepcopy(context["validation_grating"])
            replacements = []
            for path, saved in context.get("validation_profiles", {}).items():
                names = path.split(".")[2::2]
                if names == context.get("validation_path", []):
                    continue
                sibling = profile_from_dict(saved)
                sibling_compiled = compile_profile(sibling)
                target = containing
                for name in names:
                    target = target["children"][name]
                replacements.append((names, copy.deepcopy(target), sibling_compiled.node))
                input_files.update(sibling_compiled.files)
                profiles.append((sibling, sibling_compiled.fingerprint))
            target = containing
            for name in context.get("validation_path", []):
                target = target["children"][name]
            replacements.append((context.get("validation_path", []), copy.deepcopy(target), compiled.node))
            containing = remap_boundary_media(containing, [(old, new) for names, old, new in replacements])
            for names, old, new in replacements:
                target = containing
                for name in names:
                    target = target["children"][name]
                target.clear()
                target.update(copy.deepcopy(new))
            input_files.update(referenced_files(containing))
            pending = [containing]
            while pending:
                current = pending.pop()
                for name in ("medium_i", "medium_t"):
                    if name in current.get("params", {}):
                        current["params"][name] = native_material_value(current["params"][name])
                pending.extend(current.get("children", {}).values())
        signatures = {}
        for path in input_files:
            value = Path(path).stat()
            signatures[path] = (value.st_mtime_ns, value.st_size, value.st_ino)
        metadata = {"node": compiled.node, "materials": compiled.materials,
                    "dimension": document.dimension, "context": payload["context"],
                    "fingerprint": compiled.fingerprint, "file_signatures": signatures,
                    "one_d_mode": getattr(compiled, "one_d_mode", "surface"),
                    "cross_section": getattr(compiled, "cross_section", []),
                    "boundaries_1d": getattr(compiled, "boundaries_1d", [])}
        if metadata["one_d_mode"] != "native":
            np.savez_compressed(folder / "design.npz", metadata=np.array(json.dumps(metadata)),
                                faces=compiled.faces, material_ids=compiled.material_ids,
                                bounds=compiled.bounds, lattice=compiled.lattice)
            print(json.dumps({"event": "design"}), flush=True)
        if sys.stdin.readline().strip() != "validate":
            return 0

        def progress(message):
            print(json.dumps({"event": "status", "message": str(message)}), flush=True)

        request = GeometryRequest(grating=compiled.node, wavelength_um=context["wavelength_um"],
                                  incident_deg=context.get("incident_deg", 0),
                                  rotation_deg=context.get("rotation_deg", 0),
                                  configuration=context.get("configuration", 0),
                                  source_label="Custom profile")
        scene = build_scene(request, progress=progress)
        from scatmech_gratings import build_default_cross_grating

        if metadata["one_d_mode"] == "native":
            from geometry_mesh import build_mesh

            mesh = build_mesh(scene, (1, 1))
            metadata["materials"] = [(label, "") for label, tensor in mesh.materials]
            material_tensors = [tensor for label, tensor in mesh.materials]
            np.savez_compressed(folder / "design.npz", metadata=np.array(json.dumps(metadata)),
                                faces=mesh.faces, material_ids=mesh.material_ids,
                                bounds=mesh.bounds, lattice=scene.lattice)
            print(json.dumps({"event": "design"}), flush=True)
        else:
            material_tensors = [scene.medium_i, scene.medium_t]
            evaluated = {}
            for label, material in compiled.materials[2:]:
                value = native_material_value(material)
                if value not in evaluated:
                    progress("Evaluating preview material: " + str(label))
                    interface = build_default_cross_grating("Null_CrossGrating")
                    interface["params"]["medium_t"] = value
                    sample = build_scene(GeometryRequest(interface, request.wavelength_um))
                    evaluated[value] = sample.medium_t
                    scene.files = list(dict.fromkeys(scene.files + sample.files))
                material_tensors.append(evaluated[value])
        np.save(folder / "material_tensors.npy", np.asarray(material_tensors, dtype=complex), allow_pickle=False)
        if containing is not None:
            progress("Checking compatibility with the containing grating…")
            parent_request = GeometryRequest(grating=containing, wavelength_um=request.wavelength_um,
                incident_deg=request.incident_deg, rotation_deg=request.rotation_deg,
                configuration=request.configuration, source_label="Containing grating")
            try:
                parent_scene = build_scene(parent_request, progress=progress)
            except Exception as error:
                raise ValueError("The containing grating is incompatible with this draft: " + str(error)) from error
            scene.files = list(dict.fromkeys(scene.files + parent_scene.files))
        if any(compile_profile(profile).fingerprint != fingerprint for profile, fingerprint in profiles):
            raise ValueError("A material file changed during validation. Refresh to evaluate the latest contents.")
        save_scene(scene, folder / "scene.npz")
        print(json.dumps({"event": "validated"}), flush=True)
        return 0
    except Exception as error:
        print(f"{type(error).__name__}: {error}", file=sys.stderr, flush=True)
        return 1


if __name__ == "__main__":
    sys.exit(main())
