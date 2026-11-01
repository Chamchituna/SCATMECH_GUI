from __future__ import annotations

from dataclasses import asdict
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import traceback


def run(report_path):
    report_path = Path(report_path).resolve()
    report = {"passed": False, "checks": []}
    gui = None
    for name in ("stdout", "stderr"):
        if getattr(sys, name) is None:
            setattr(sys, name, open(os.devnull, "w", encoding="utf-8"))
    try:
        from scatmech_paths import configure_runtime, get_data_dir, worker_command
        configure_runtime(change_directory=True)
        from main import QApplication, SCATMECHGui
        from custom_profile import default_editor_profile, profile_to_dict
        from geometry_scene import load_scene
        from geometry_types import GeometryRequest
        from scatmech_gratings import build_default_one_d_grating

        app = QApplication.instance() or QApplication([])
        gui = SCATMECHGui()
        report["checks"].append("All four GUI pages constructed")
        gui.brdf_form.azimuth_step.setText("30")
        gui.mie_form.step_angle_deg.setText("30")
        gui.mie_form.step_phi_deg.setText("30")
        gui.rcw_form.order.setText("2")
        gui.rcw_form.sweep_start_um.setText("0.532")
        gui.rcw_form.sweep_end_um.setText("0.542")
        gui.rcw_form.sweep_step_um.setText("0.01")
        for name, form, method in (
            ("BRDF", gui.brdf_form, "run_brdfprog"),
            ("Reflect", gui.reflect_form, "run_reflectprog"),
            ("Mie", gui.mie_form, "run_mieprog"),
            ("RCW", gui.rcw_form, "run_rcwprog"),
        ):
            getattr(form, method)()
            output = getattr(form, "last_csv_path", None)
            if not output or not Path(output).is_file() or not getattr(form, "last_plot_success", False):
                widget = getattr(form, "log", None) or getattr(form, "output_box", None)
                detail = widget.toPlainText() if widget is not None else "No output was produced."
                raise RuntimeError(f"{name} solver/plot failed: {detail}")
            form.figure.savefig(get_data_dir() / f"smoke-{name.lower()}.png")
            report["checks"].append(f"{name}: native calculation, CSV parsing and Matplotlib rendering")
        app.processEvents()


        with tempfile.TemporaryDirectory(prefix="SCATMECH smoke Ω ") as temporary:
            folder = Path(temporary)
            request_path = folder / "geometry request.json"
            scene_path = folder / "geometry scene.npz"
            request = GeometryRequest(build_default_one_d_grating("Single_Line_Grating"), .532, resolution=16)
            request_path.write_text(json.dumps(asdict(request)), encoding="utf-8")
            program, arguments = worker_command(Path(__file__).with_name("geometry_worker.py"), request_path, scene_path)
            result = subprocess.run([program, *arguments], input="", capture_output=True, text=True,
                                    timeout=60, creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0))
            if result.returncode or not scene_path.is_file():
                raise RuntimeError(f"Geometry worker failed: {result.stdout}\n{result.stderr}")
            scene = load_scene(scene_path)
            if not scene.layers:
                raise RuntimeError("Geometry worker returned an empty scene.")
            report["checks"].append("Geometry worker: native helper, JSON progress, saved geometry; spaces and Unicode paths")

            for dimension in (1, 2):
                output_folder = folder / f"custom profile {dimension}"
                output_folder.mkdir()
                profile_request = output_folder / "request.json"
                profile_request.write_text(json.dumps({
                    "document": profile_to_dict(default_editor_profile(dimension)),
                    "context": {"wavelength_um": .532},
                }), encoding="utf-8")
                program, arguments = worker_command(Path(__file__).with_name("custom_profile_worker.py"), profile_request, output_folder)
                result = subprocess.run([program, *arguments], input="validate\n", capture_output=True,
                                        text=True, timeout=60,
                                        creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0))
                events = [json.loads(line) for line in result.stdout.splitlines() if line.startswith("{")]
                if result.returncode or not any(event.get("event") == "validated" for event in events):
                    raise RuntimeError(f"Custom profile worker {dimension}D failed: {result.stdout}\n{result.stderr}")
                load_scene(output_folder / "scene.npz")
                report["checks"].append(f"Custom profile {dimension}D: stdin validation, native geometry and saved scene")
        report["passed"] = True
        return 0
    except Exception:
        report["error"] = traceback.format_exc()
        return 1
    finally:
        if gui is not None:
            gui.close()
        report_path.parent.mkdir(parents=True, exist_ok=True)
        report_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
