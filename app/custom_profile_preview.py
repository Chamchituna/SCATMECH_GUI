from __future__ import annotations

import copy
import json
import os
from pathlib import Path
import shutil
import signal
import sys
import tempfile
from scatmech_paths import worker_command
from types import SimpleNamespace

import numpy as np
from PyQt5.QtCore import QProcess, Qt, QTimer, pyqtSignal
from PyQt5.QtWidgets import QApplication, QComboBox, QHBoxLayout, QLabel, QPushButton, QSpinBox, QToolButton, QVBoxLayout, QWidget
from matplotlib import colormaps, patheffects
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib.colors import LightSource
from matplotlib.figure import Figure
from matplotlib.patches import Patch
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

from geometry_mesh import build_mesh
from geometry_navigation import ConstrainedAxes3D
from geometry_scene import load_scene
from geometry_viewer import incident_direction



def preview_resolution_summary(document, scene, mode, face_count, bounds):
    if mode == 1 and scene is not None:
        if scene.dimension == 2:
            nu = max((len(layer.u_edges) - 1 for layer in scene.layers), default=0)
            nv = max((len(layer.v_edges) - 1 for layer in scene.layers), default=0)
            text = f"Native grid {document.grid1}×{document.grid2}; preview up to {nu}×{nv}; {len(scene.layers)} layers retained. Sampled preview; narrow features may be omitted."
        else:
            text = f"{len(scene.layers)} native layers retained · {face_count:,} displayed faces."
    else:
        text = f"Design boundaries · {face_count:,} faces."
        if document is not None and document.dimension == 2 and any(feature.get("shape") == "circle" for feature in document.features):
            text += f" Circles: {document.curve_segments} segments."
    if document is not None and document.dimension == 1:
        text += f" Invariant y: {np.ptp(bounds[1]):.6g} µm (finite illustration)."
    return text


class CustomProfilePreview(QWidget):
    validation_changed = pyqtSignal(int, bool, str)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.document = None
        self.compiled = None
        self.revision = 0
        self.validated_revision = -1
        self.worker_path = Path(__file__).with_name("custom_profile_worker.py")
        self._context = {}
        self._jobs = {}
        self._scene = None
        self._material_tensors = None
        self._design_revision = -1
        self._bounds = None
        self._camera_initialized = False
        self._closed = False
        self._file_signatures = {}
        self.file_poll = QTimer(self)
        self.file_poll.setInterval(1000)
        self.file_poll.timeout.connect(self._poll_files)
        self.debounce = QTimer(self)
        self.debounce.setSingleShot(True)
        self.debounce.setInterval(500)
        self.debounce.timeout.connect(self._validate)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        controls = QHBoxLayout()
        self.mode = QComboBox()
        self.mode.addItems(["Design", "Solver layers"])
        self.mode.currentIndexChanged.connect(self._render)
        controls.addWidget(self.mode)
        self.refresh_button = QPushButton("Refresh")
        self.refresh_button.clicked.connect(self.refresh)
        controls.addWidget(self.refresh_button)
        self.cancel_button = QPushButton("Cancel update")
        self.cancel_button.clicked.connect(self.cancel)
        self.cancel_button.setEnabled(False)
        controls.addWidget(self.cancel_button)
        controls.addStretch()
        layout.addLayout(controls)
        views = QHBoxLayout()
        for label, callback in (("Fit", self.fit), ("Reset", self.reset_view),
                                ("Top", lambda: self.set_view(90, -90)),
                                ("Front", lambda: self.set_view(0, -90)),
                                ("Side", lambda: self.set_view(0, 0))):
            button = QPushButton(label)
            button.setFixedWidth(48)
            button.clicked.connect(callback)
            views.addWidget(button)
        views.addWidget(QLabel("Repeats"))
        self.repeat_x = QSpinBox()
        self.repeat_y = QSpinBox()
        for widget in (self.repeat_x, self.repeat_y):
            widget.setRange(1, 5)
            widget.setFixedWidth(55)
            widget.setValue(3)
            widget.valueChanged.connect(self._render)
            views.addWidget(widget)
        views.addStretch()
        layout.addLayout(views)
        self.figure = Figure(figsize=(6, 5))
        self.canvas = FigureCanvasQTAgg(self.figure)
        self.axes = self.figure.add_subplot(111, axes_class=ConstrainedAxes3D, computed_zorder=False)
        self.axes.geometry_bounds = lambda: self._bounds
        self.figure.subplots_adjust(left=.08, right=.92, top=.87, bottom=.12)
        self.axes.set_anchor("C")
        self.canvas.mpl_connect("scroll_event", self._scroll)
        layout.addWidget(NavigationToolbar2QT(self.canvas, self))
        layout.addWidget(self.canvas, 1)
        self.resolution_label = QLabel("Design preview")
        self.resolution_label.setWordWrap(True)
        self.resolution_label.setTextFormat(Qt.PlainText)
        resolution_row = QHBoxLayout()
        resolution_row.addWidget(self.resolution_label, 1)
        self.details_button = QToolButton()
        self.details_button.setText("Preview details")
        self.details_button.setCheckable(True)
        self.details_button.setArrowType(Qt.RightArrow)
        self.details_button.setToolButtonStyle(Qt.ToolButtonTextBesideIcon)
        self.details_button.toggled.connect(self._toggle_details)
        resolution_row.addWidget(self.details_button, 0, Qt.AlignTop)
        layout.addLayout(resolution_row)
        self.details_label = QLabel()
        self.details_label.setWordWrap(True)
        self.details_label.setTextFormat(Qt.PlainText)
        self.details_label.setTextInteractionFlags(Qt.TextSelectableByMouse)
        self.details_label.hide()
        layout.addWidget(self.details_label)
        self.status_label = QLabel("Ready")
        self.status_label.setWordWrap(True)
        self.status_label.setTextFormat(Qt.PlainText)
        layout.addWidget(self.status_label)
        self._stale = self.figure.text(.5, .96, "", ha="center", va="top", color="#a52c22", fontsize=10)
        self._title = self.figure.text(.5, .88, "", ha="center", va="top", fontsize=11)
        app = QApplication.instance()
        if app is not None:
            app.aboutToQuit.connect(self.shutdown)

    def _toggle_details(self, expanded):
        self.details_label.setVisible(expanded)
        self.details_button.setArrowType(Qt.DownArrow if expanded else Qt.RightArrow)

    def _status(self, message, valid=False):
        self.status_label.setText(message)
        self.validation_changed.emit(self.revision, valid, message)

    def set_document(self, document, context):
        from custom_profile import profile_to_dict

        self._closed = False
        self.file_poll.start()
        self.document = copy.deepcopy(document)
        self._context = copy.deepcopy(context)
        native = document.dimension == 1 and getattr(document, "one_d_mode", "surface") == "native"
        blocked = self.mode.blockSignals(True)
        self.mode.setItemText(0, "Native geometry" if native else "Design")
        self.mode.setEnabled(not native)
        if native:
            self.mode.setCurrentIndex(0)
        self.mode.blockSignals(blocked)
        self.revision += 1
        self.validated_revision = -1
        self.debounce.stop()
        self._cancel_jobs()
        self._stale.set_text("Updating design…")
        self._status("Updating design; validation pending.")
        self.canvas.draw_idle()
        directory = None
        try:
            payload = {"document": profile_to_dict(document), "context": self._context}
            serialized = json.dumps(payload, allow_nan=False)
            directory = Path(tempfile.mkdtemp(prefix="scatmech-profile-preview-"))
            request = directory / "request.json"
            request.write_text(serialized, encoding="utf-8")
            process = QProcess(self)
            job = {"revision": self.revision, "directory": directory, "process": process,
                   "pid": 0, "group_owned": False, "cancelled": False,
                   "stdout": bytearray(), "stderr": bytearray(), "validate": False}
            self._jobs[self.revision] = job
            process.started.connect(lambda: self._started(job))
            process.readyReadStandardOutput.connect(lambda: self._stdout(job))
            process.readyReadStandardError.connect(lambda: self._stderr(job))
            process.finished.connect(lambda code, status: self._finished(job, code, status))
            process.errorOccurred.connect(lambda error: self._error(job, error))
            program, arguments = worker_command(self.worker_path, request, directory)
            process.start(program, arguments)
            self.cancel_button.setEnabled(True)
            self.debounce.start()
        except Exception as error:
            if directory is not None:
                shutil.rmtree(directory, ignore_errors=True)
            self._stale.set_text("Invalid draft — previous preview")
            self._status(str(error))

    def _started(self, job):
        job["pid"] = int(job["process"].processId())
        self._check_group(job)
        if job["cancelled"]:
            self._stop(job)
        elif job["validate"]:
            job["process"].write(b"validate\n")

    def _check_group(self, job):
        if os.name == "posix" and job["pid"]:
            try:
                job["group_owned"] = os.getpgid(job["pid"]) == job["pid"]
            except ProcessLookupError:
                pass

    def _validate(self):
        job = self._jobs.get(self.revision)
        if job is not None and not job["cancelled"]:
            job["validate"] = True
            if job["process"].state() == QProcess.Running:
                job["process"].write(b"validate\n")

    def _current(self, job):
        return job["revision"] == self.revision and not job["cancelled"] and not self._closed

    def _stdout(self, job):
        job["stdout"].extend(bytes(job["process"].readAllStandardOutput()))
        self._check_group(job)
        while b"\n" in job["stdout"]:
            line, job["stdout"] = job["stdout"].split(b"\n", 1)
            if not self._current(job):
                continue
            try:
                event = json.loads(line)
                if event.get("event") == "design":
                    with np.load(job["directory"] / "design.npz", allow_pickle=False) as archive:
                        metadata = json.loads(archive["metadata"].item())
                        self.compiled = SimpleNamespace(**metadata, **{key: archive[key].copy() for key in
                            ("faces", "material_ids", "bounds", "lattice")})
                    self._design_revision = self.revision
                    self._stale.set_text("Native geometry · solver validation pending" if getattr(self.compiled, "one_d_mode", "surface") == "native"
                                         else "Design ready · solver validation pending")
                    self._render()
                elif event.get("event") == "status":
                    self.status_label.setText(str(event.get("message", "Validating…")))
            except Exception as error:
                self._status("Preview failed: " + str(error))

    def _stderr(self, job):
        job["stderr"].extend(bytes(job["process"].readAllStandardError()))
        del job["stderr"][:-64000]

    def _error(self, job, error):
        if error == QProcess.FailedToStart:
            job["stderr"].extend(job["process"].errorString().encode())
            self._finished(job, -1, QProcess.CrashExit)

    def _finished(self, job, code, status):
        if job["revision"] not in self._jobs:
            return
        self._stdout(job)
        self._stderr(job)
        if self._current(job):
            try:
                if code != 0 or status != QProcess.NormalExit:
                    raise RuntimeError(job["stderr"].decode(errors="replace").strip() or "Native validation failed.")
                if self._design_revision != self.revision:
                    raise RuntimeError("The current design preview was not generated.")
                self._scene = load_scene(job["directory"] / "scene.npz")
                material_file = job["directory"] / "material_tensors.npy"
                self._material_tensors = np.load(material_file, allow_pickle=False) if material_file.is_file() else None
                self._file_signatures = {path: self._file_signature(path) for path in self._scene.files}
                originals = getattr(self.compiled, "file_signatures", {})
                for path, signature in originals.items():
                    if self._file_signature(path) != tuple(signature):
                        raise ValueError("A material file changed during validation. Refresh to evaluate the latest contents.")
                    self._file_signatures[path] = tuple(signature)
                self.validated_revision = self.revision
                self._stale.set_text("")
                if not self._render():
                    raise RuntimeError(self.status_label.text())
                self._status("Validated — ready to apply to simulation.", True)
            except Exception as error:
                self.validated_revision = -1
                self._stale.set_text("Draft not validated")
                self._status(str(error)[-3000:])
                self.canvas.draw_idle()
        if job["cancelled"]:
            self._force_kill(job)
        self._jobs.pop(job["revision"], None)
        shutil.rmtree(job["directory"], ignore_errors=True)
        job["process"].deleteLater()
        self.cancel_button.setEnabled(any(not item["cancelled"] for item in self._jobs.values()))

    def refresh(self):
        if self.document is not None:
            self.set_document(self.document, self._context)

    @staticmethod
    def _file_signature(path):
        try:
            value = Path(path).stat()
            return value.st_mtime_ns, value.st_size, value.st_ino
        except OSError:
            return None

    def is_validated(self):
        return self.validated_revision == self.revision and all(
            self._file_signature(path) == signature for path, signature in self._file_signatures.items())

    def _poll_files(self):
        if self.validated_revision == self.revision and not self.is_validated():
            self.refresh()

    def cancel(self):
        self.debounce.stop()
        self.file_poll.stop()
        self.revision += 1
        self.validated_revision = -1
        self._cancel_jobs()
        self._stale.set_text("Update cancelled")
        self._status("Update cancelled. Refresh to validate this draft.")
        self.canvas.draw_idle()

    def _cancel_jobs(self):
        for job in list(self._jobs.values()):
            if not job["cancelled"]:
                job["cancelled"] = True
                self._stop(job)
        self.cancel_button.setEnabled(False)

    def _stop(self, job):
        self._check_group(job)
        if job["group_owned"]:
            try:
                os.killpg(job["pid"], signal.SIGTERM)
            except ProcessLookupError:
                pass
        elif job["process"].state() != QProcess.NotRunning:
            job["process"].terminate()
        QTimer.singleShot(1000, lambda: self._force_kill(job) if job["revision"] in self._jobs else None)

    def _force_kill(self, job):
        if job["group_owned"]:
            try:
                os.killpg(job["pid"], signal.SIGKILL)
            except ProcessLookupError:
                pass
        if job["process"].state() != QProcess.NotRunning:
            job["process"].kill()

    def shutdown(self):
        self._closed = True
        self.canvas._draw_pending = False
        self.debounce.stop()
        self.file_poll.stop()
        self._cancel_jobs()
        for job in list(self._jobs.values()):
            if job["process"].state() != QProcess.NotRunning and not job["process"].waitForFinished(1000):
                self._force_kill(job)
                job["process"].waitForFinished(500)

    def _camera(self):
        return (self.axes.elev, self.axes.azim, getattr(self.axes, "roll", 0),
                self.axes.get_xlim3d(), self.axes.get_ylim3d(), self.axes.get_zlim3d())

    def _design_mesh(self):
        from custom_profile_mesh import build_design_mesh

        return build_design_mesh(self.compiled, (self.repeat_x.value(), self.repeat_y.value()))

    def _render(self, *_):
        try:
            self._render_scene()
            return True
        except Exception as error:
            self.validated_revision = -1
            self._stale.set_text("Preview could not be rendered")
            self._status("Preview failed: " + str(error))
            return False

    def _render_scene(self):
        if self.compiled is None:
            return
        camera = self._camera()
        native = getattr(self.compiled, "one_d_mode", "surface") == "native"
        solver = self.mode.currentIndex() == 1 or native
        if solver:
            if self._scene is None or self.validated_revision != self.revision:
                self._stale.set_text("Solver layers out of date · showing previous preview")
                self.canvas.draw_idle()
                return
            mesh = build_mesh(self._scene, (self.repeat_x.value(), self.repeat_y.value()))
            faces, ids, materials, bounds = mesh.faces, mesh.material_ids, mesh.materials, mesh.bounds
            if self._material_tensors is not None:
                mapping = []
                materials = list(self.compiled.materials)
                for label, tensor in mesh.materials:
                    matches = np.flatnonzero(np.all(np.isclose(self._material_tensors, tensor, rtol=1e-9, atol=1e-12), axis=1))
                    if len(matches):
                        mapping.append(int(matches[0]))
                    else:
                        mapping.append(len(materials))
                        materials.append((label, tensor))
                ids = np.asarray(mapping, dtype=int)[ids]
            notes = " · ".join(mesh.notes)
        else:
            faces, ids, materials, bounds = self._design_mesh()
            notes = "Design boundaries; native layer/grid approximation is available under Solver layers."
            if self.compiled.dimension == 2 and self.document is not None:
                notes += f" Circles use {self.document.curve_segments} boundary segments."
        visible = ids != 0
        self.axes.clear()
        self.axes.set_anchor("C")
        colors = [colormaps["tab10"](((int(index) - 1) % 10) / 9) for index in ids[visible]]
        if np.any(visible):
            surface = Poly3DCollection(faces[visible], facecolors=colors, linewidths=0,
                                       antialiased=False, shade=True, lightsource=LightSource(315, 45), zorder=2)
            surface.set_edgecolors("none")
            self.axes.add_collection3d(surface)
        self._bounds = bounds
        self.axes.set_xlabel("x (µm)")
        self.axes.set_ylabel("y (µm)")
        self.axes.set_zlabel("z (µm)")
        self._title.set_text(f"Custom {self.compiled.dimension}D profile · {self.mode.currentText()}")
        direction = incident_direction(SimpleNamespace(incident_deg=self._context.get("incident_deg", 0),
                    rotation_deg=self._context.get("rotation_deg", 0), configuration=self._context.get("configuration", 0)))
        length = max(float(np.max(np.ptp(bounds, axis=1))) * .25, .05)
        target = np.array([np.mean(bounds[0]), np.mean(bounds[1]), bounds[2, 1]])
        if self._context.get("configuration", 0) >= 2:
            target[2] = bounds[2, 0]
        self._beam_start = target - length * direction
        beam = self.axes.quiver(*self._beam_start, *direction, length=length, arrow_length_ratio=.22,
                               color="#c63b26", linewidth=2.5, zorder=20)
        beam.set_path_effects([patheffects.Stroke(linewidth=4, foreground="white"), patheffects.Normal()])
        self.axes.text(*self._beam_start, " Incident beam", color="#b42f20", fontsize=9, zorder=21,
                       bbox={"facecolor": "white", "edgecolor": "none", "alpha": .85, "pad": 2})
        handles = [Patch(facecolor=colormaps["tab10"](((index - 1) % 10) / 9), label=str(material[0]))
                   for index, material in enumerate(materials) if index and np.any(ids[visible] == index)]
        if handles:
            self.axes.legend(handles=handles, loc="upper left", fontsize=8,
                             bbox_to_anchor=(.02, .81), bbox_transform=self.figure.transFigure)
        if self._camera_initialized:
            self.axes.view_init(elev=camera[0], azim=camera[1], roll=camera[2])
            self.axes.set_xlim3d(camera[3])
            self.axes.set_ylim3d(camera[4])
            self.axes.set_zlim3d(camera[5])
            self._aspect()
        else:
            self.axes.view_init(elev=28, azim=-58)
            self.fit()
            self._camera_initialized = True
        if not solver and self._design_revision == self.revision:
            self._stale.set_text("" if self.validated_revision == self.revision else "Design ready · solver validation pending")
        text = f"{len(faces):,} displayed faces. " + notes
        if self.compiled.dimension == 1:
            text += f" Invariant y direction shown over a finite {np.ptp(bounds[1]):.6g} µm length."
        self.details_label.setText(text)
        self.details_button.setToolTip(text)
        self.resolution_label.setToolTip(text)
        self.resolution_label.setText(preview_resolution_summary(
            self.document, self._scene, int(solver), len(faces), bounds))
        self.canvas.draw_idle()

    def _aspect(self):
        extent = np.array([np.ptp(self.axes.get_xlim3d()), np.ptp(self.axes.get_ylim3d()), np.ptp(self.axes.get_zlim3d())])
        self.axes.set_box_aspect(np.maximum(extent, max(extent) * .05))

    def fit(self, *_):
        if self._bounds is None:
            return
        bounds = self._bounds.copy()
        bounds[:, 0] = np.minimum(bounds[:, 0], self._beam_start)
        bounds[:, 1] = np.maximum(bounds[:, 1], self._beam_start)
        margin = np.maximum(np.ptp(bounds, axis=1), .01) * .08
        for index, setter in enumerate((self.axes.set_xlim3d, self.axes.set_ylim3d, self.axes.set_zlim3d)):
            setter(bounds[index, 0] - margin[index], bounds[index, 1] + margin[index])
        self._aspect()
        self.canvas.draw_idle()

    def set_view(self, elevation, azimuth):
        self.axes.view_init(elev=elevation, azim=azimuth)
        self.canvas.draw_idle()

    def reset_view(self):
        self.set_view(28, -58)
        self.fit()

    def _scroll(self, event):
        if event.inaxes != self.axes:
            return
        scale = .85 if event.button == "up" else 1 / .85
        for getter, setter in ((self.axes.get_xlim3d, self.axes.set_xlim3d), (self.axes.get_ylim3d, self.axes.set_ylim3d),
                               (self.axes.get_zlim3d, self.axes.set_zlim3d)):
            low, high = getter()
            middle, half = (low + high) / 2, (high - low) * scale / 2
            setter(middle - half, middle + half)
        self.axes.constrain_pan()
        self._aspect()
        self.canvas.draw_idle()
