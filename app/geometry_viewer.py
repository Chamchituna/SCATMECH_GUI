from __future__ import annotations

from dataclasses import asdict, replace
import json
import math
import os
from pathlib import Path
import shutil
import signal
import sys
import tempfile

from scatmech_paths import worker_command

import numpy as np
from PyQt5.QtCore import QFileSystemWatcher, QProcess, Qt, QTimer
from PyQt5.QtGui import QColor
from PyQt5.QtWidgets import (
    QApplication, QCheckBox, QComboBox, QDoubleSpinBox, QFileDialog,
    QHBoxLayout, QHeaderView, QLabel, QPushButton, QSlider, QSpinBox,
    QSplitter, QTableWidget, QTableWidgetItem, QVBoxLayout, QWidget,
)
from matplotlib import colormaps
from matplotlib import patheffects
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib.colors import LightSource, Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.figure import Figure
from matplotlib.patches import Patch
from matplotlib.ticker import AutoLocator, NullLocator
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

from geometry_mesh import build_mesh
from geometry_navigation import ConstrainedAxes3D
from geometry_scene import load_scene
from geometry_types import GeometryRequest


def incident_direction(scene):
    theta = math.radians(scene.incident_deg)
    phi = math.radians(scene.rotation_deg)
    return np.array([
        math.sin(theta) * math.cos(phi),
        -math.sin(theta) * math.sin(phi),
        math.cos(theta) * (1 if scene.configuration >= 2 else -1),
    ])


def request_files(request):
    from geometry_native import referenced_files

    found = set(referenced_files(request.grating))

    def visit(value):
        if isinstance(value, dict):
            for item in value.values():
                visit(item)
        elif isinstance(value, (list, tuple)):
            for item in value:
                visit(item)
        elif isinstance(value, str):
            candidate = Path(value.strip().strip('"').partition("!")[0]).expanduser()
            try:
                if candidate.is_file() or (candidate.is_absolute() and candidate.parent.is_dir()):
                    found.add(str(candidate.resolve()))
            except (OSError, ValueError):
                pass

    visit(request.grating)
    return found


class GeometryViewer(QWidget):
    def __init__(self, parent=None, request_provider=None):
        super().__init__(parent, Qt.Window)
        self.request_provider = request_provider
        self.worker_path = Path(__file__).with_name("geometry_worker.py")
        self.scene = None
        self.mesh = None
        self.stale = True
        self._generation = 0
        self._jobs = {}
        self._closed = False
        self._watched_files = set()
        self._file_signatures = {}
        self._material_state = {}
        self._table_keys = []
        self._color_selected = False
        self._camera_initialized = False
        self._colorbar = None
        self._last_request = None
        self.setWindowTitle("3D RCW Geometry")
        self.resize(1150, 900)
        self.setMinimumSize(680, 520)
        self.debounce = QTimer(self)
        self.debounce.setSingleShot(True)
        self.debounce.setInterval(500)
        self.debounce.timeout.connect(self.refresh)
        self.watcher = QFileSystemWatcher(self)
        self.watcher.fileChanged.connect(self._file_changed)
        self.watcher.directoryChanged.connect(self._directory_changed)
        self.file_poll = QTimer(self)
        self.file_poll.setInterval(1000)
        self.file_poll.timeout.connect(self._poll_files)
        self._build_ui()
        application = QApplication.instance()
        if application is not None:
            application.aboutToQuit.connect(self.shutdown)

    def _build_ui(self):
        layout = QVBoxLayout(self)
        actions = QHBoxLayout()
        self.refresh_button = QPushButton("Refresh")
        self.refresh_button.clicked.connect(self.refresh)
        actions.addWidget(self.refresh_button)
        self.cancel_button = QPushButton("Cancel")
        self.cancel_button.clicked.connect(self.cancel)
        self.cancel_button.setEnabled(False)
        actions.addWidget(self.cancel_button)
        self.wavelength_override = QCheckBox("Preview wavelength (µm)")
        self.wavelength_override.toggled.connect(self._override_changed)
        actions.addWidget(self.wavelength_override)
        self.wavelength = QDoubleSpinBox()
        self.wavelength.setDecimals(8)
        self.wavelength.setRange(0.00000001, 1000000)
        self.wavelength.setValue(0.532)
        self.wavelength.setEnabled(False)
        self.wavelength.valueChanged.connect(self.source_changed)
        actions.addWidget(self.wavelength)
        actions.addStretch()
        self.export_button = QPushButton("Export PNG")
        self.export_button.clicked.connect(lambda: self.export_png())
        self.export_button.setEnabled(False)
        actions.addWidget(self.export_button)
        layout.addLayout(actions)
        view_controls = QHBoxLayout()
        for label, callback in (
            ("Fit", self.fit), ("Reset", self.reset_view),
            ("Top", lambda: self.set_view(90, -90)),
            ("Front", lambda: self.set_view(0, -90)),
            ("Side", lambda: self.set_view(0, 0)),
        ):
            button = QPushButton(label)
            button.clicked.connect(callback)
            view_controls.addWidget(button)
        view_controls.addWidget(QLabel("Repeat a₁"))
        self.repeat_x = QSpinBox()
        self.repeat_x.setRange(1, 5)
        self.repeat_x.setValue(3)
        self.repeat_x.valueChanged.connect(self._rebuild_mesh)
        view_controls.addWidget(self.repeat_x)
        self.repeat_y_label = QLabel("Repeat a₂")
        view_controls.addWidget(self.repeat_y_label)
        self.repeat_y = QSpinBox()
        self.repeat_y.setRange(1, 5)
        self.repeat_y.setValue(3)
        self.repeat_y.valueChanged.connect(self._rebuild_mesh)
        view_controls.addWidget(self.repeat_y)
        view_controls.addStretch()
        layout.addLayout(view_controls)
        coloring = QHBoxLayout()
        coloring.addWidget(QLabel("Color"))
        self.color_mode = QComboBox()
        self.color_mode.addItem("Material", None)
        for component, symbol in enumerate(("εₓ", "εᵧ", "εz", "μₓ", "μᵧ", "μz")):
            self.color_mode.addItem("Re " + symbol, (component, "real"))
            self.color_mode.addItem("Im " + symbol, (component, "imag"))
        self.color_mode.currentIndexChanged.connect(self._color_changed)
        coloring.addWidget(self.color_mode)
        self.smooth_cylinders = QCheckBox("Smooth cylinders")
        self.smooth_cylinders.setChecked(True)
        self.smooth_cylinders.setToolTip("Draw the native radius at each solver layer. Turn off to inspect sampled grid cells.")
        self.smooth_cylinders.toggled.connect(self._rebuild_mesh)
        self.smooth_cylinders.hide()
        coloring.addWidget(self.smooth_cylinders)
        coloring.addWidget(QLabel("Cut away above"))
        self.slice_axis = QComboBox()
        self.slice_axis.addItems(["None", "X", "Y", "Z"])
        self.slice_axis.currentIndexChanged.connect(self._slice_changed)
        coloring.addWidget(self.slice_axis)
        self.slice_slider = QSlider(Qt.Horizontal)
        self.slice_slider.setRange(1, 100)
        self.slice_slider.setValue(100)
        self.slice_slider.setEnabled(False)
        self.slice_slider.setTracking(False)
        self.slice_slider.valueChanged.connect(self._slice_changed)
        coloring.addWidget(self.slice_slider, 1)
        self.slice_label = QLabel("100% retained")
        coloring.addWidget(self.slice_label)
        layout.addLayout(coloring)
        splitter = QSplitter(Qt.Vertical)
        plot = QWidget()
        plot_layout = QVBoxLayout(plot)
        plot_layout.setContentsMargins(0, 0, 0, 0)
        self.figure = Figure(figsize=(9, 6))
        self.canvas = FigureCanvasQTAgg(self.figure)
        self.axes = self.figure.add_subplot(111, axes_class=ConstrainedAxes3D, computed_zorder=False)
        self.axes.geometry_bounds = lambda: self.mesh.bounds if self.mesh is not None else None
        self.figure.subplots_adjust(left=0.10, right=0.90, top=0.89, bottom=0.08)
        self.axes.set_anchor("C")
        self.stale_artist = self.figure.text(
            0.5, 0.97, "OUT OF DATE", ha="center", va="top", color="#b53030",
            fontsize=12, weight="bold", bbox={"facecolor": "white", "alpha": 0.9, "edgecolor": "none"},
        )
        self.toolbar = NavigationToolbar2QT(self.canvas, self)
        self.toolbar.setToolTip("Drag to rotate; use Pan and Zoom or the mouse wheel to change the view.")
        self.canvas.mpl_connect("scroll_event", self._scroll)
        self.canvas.mpl_connect("button_release_event", self._camera_changed)
        plot_layout.addWidget(self.toolbar)
        plot_layout.addWidget(self.canvas, 1)
        splitter.addWidget(plot)
        self.material_table = QTableWidget(0, 3)
        self.material_table.setHorizontalHeaderLabels(["Material / visibility", "Opacity", "Complex ε and μ components"])
        self.material_table.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeToContents)
        self.material_table.horizontalHeader().setSectionResizeMode(1, QHeaderView.ResizeToContents)
        self.material_table.horizontalHeader().setSectionResizeMode(2, QHeaderView.Stretch)
        self.material_table.verticalHeader().setVisible(False)
        self.material_table.itemChanged.connect(self._material_changed)
        splitter.addWidget(self.material_table)
        splitter.setStretchFactor(0, 5)
        splitter.setStretchFactor(1, 1)
        splitter.setSizes([570, 140])
        layout.addWidget(splitter, 1)
        self.resolution_label = QLabel("Solver geometry will appear here.")
        self.resolution_label.setWordWrap(True)
        self.resolution_label.setTextFormat(Qt.PlainText)
        layout.addWidget(self.resolution_label)
        self.status_label = QLabel("Ready")
        self.status_label.setTextFormat(Qt.PlainText)
        self.status_label.setWordWrap(True)
        self.status_label.setTextInteractionFlags(Qt.TextSelectableByMouse)
        layout.addWidget(self.status_label)

    def showEvent(self, event):
        super().showEvent(event)
        self._closed = False
        self.file_poll.start()
        self.refresh()

    def closeEvent(self, event):
        self._closed = True
        self.debounce.stop()
        self.file_poll.stop()
        self._generation += 1
        self._cancel_jobs()
        super().closeEvent(event)

    def shutdown(self):
        self._closed = True
        self.debounce.stop()
        self.file_poll.stop()
        self._cancel_jobs()
        for job in list(self._jobs.values()):
            process = job["process"]
            if process.state() != QProcess.NotRunning and not process.waitForFinished(1000):
                self._force_kill(job)
                process.waitForFinished(500)

    def source_changed(self, *_):
        self._generation += 1
        self._mark_stale("OUT OF DATE — Inputs changed; waiting to update…")
        self._cancel_jobs()
        if not self._closed and self.isVisible():
            self.debounce.start()

    def _override_changed(self, enabled):
        self.wavelength.setEnabled(enabled)
        self.source_changed()

    def _mark_stale(self, message):
        self.stale = True
        self.stale_artist.set_visible(True)
        self.status_label.setText(message)
        self.canvas.draw_idle()

    def refresh(self):
        if self._closed:
            return
        self.debounce.stop()
        self._generation += 1
        generation = self._generation
        self._cancel_jobs()
        self._mark_stale("OUT OF DATE — Reading current inputs…")
        try:
            if self.request_provider is None:
                raise ValueError("No source form is connected to this viewer.")
            request = self.request_provider()
            if not isinstance(request, GeometryRequest):
                raise ValueError("The source form did not return a geometry request.")
            if self.wavelength_override.isChecked():
                request = replace(request, wavelength_um=self.wavelength.value())
            else:
                self.wavelength.blockSignals(True)
                self.wavelength.setValue(request.wavelength_um)
                self.wavelength.blockSignals(False)
            if not math.isfinite(request.wavelength_um) or request.wavelength_um <= 0:
                raise ValueError("Preview wavelength must be positive and finite.")
            self._last_request = request
            self.setWindowTitle("{} — 3D RCW Geometry".format(request.source_label))
            self._watch(request_files(request))
            directory = Path(tempfile.mkdtemp(prefix="scatmech-geometry-"))
            input_path = directory / "request.json"
            output_path = directory / "scene.npz"
            try:
                input_path.write_text(json.dumps(asdict(request), allow_nan=False), encoding="utf-8")
            except Exception:
                shutil.rmtree(directory, ignore_errors=True)
                raise
            process = QProcess(self)
            process.setProcessChannelMode(QProcess.SeparateChannels)
            job = {
                "process": process, "directory": directory, "output": output_path,
                "cancelled": False, "pid": 0, "group_owned": False,
                "stderr": bytearray(), "stdout": bytearray(), "generation": generation,
            }
            self._jobs[generation] = job
            process.started.connect(lambda: self._started(job))
            process.readyReadStandardOutput.connect(lambda: self._read_stdout(job))
            process.readyReadStandardError.connect(lambda: self._read_stderr(job))
            process.finished.connect(lambda code, status: self._finished(job, code, status))
            process.errorOccurred.connect(lambda error: self._process_error(job, error))
            self.cancel_button.setEnabled(True)
            self.status_label.setText("OUT OF DATE — Evaluating solver geometry…")
            program, arguments = worker_command(self.worker_path, input_path, output_path)
            process.start(program, arguments)
        except Exception as error:
            self._mark_stale("OUT OF DATE — " + str(error))
            self.cancel_button.setEnabled(False)

    def _started(self, job):
        job["pid"] = int(job["process"].processId())
        self._check_group(job)
        if job["cancelled"]:
            self._stop_job(job)

    def _check_group(self, job):
        if os.name != "nt" and job["pid"]:
            try:
                job["group_owned"] = os.getpgid(job["pid"]) == job["pid"]
            except ProcessLookupError:
                pass

    def _read_stdout(self, job):
        job["stdout"].extend(bytes(job["process"].readAllStandardOutput()))
        self._check_group(job)
        while b"\n" in job["stdout"]:
            line, remainder = job["stdout"].split(b"\n", 1)
            job["stdout"] = remainder
            try:
                event = json.loads(line.decode("utf-8", errors="replace"))
            except (ValueError, TypeError):
                continue
            if job["generation"] == self._generation and not job["cancelled"]:
                if isinstance(event, dict) and event.get("event") == "status":
                    self.status_label.setText("OUT OF DATE — " + str(event.get("message", "Working…")))

    def _read_stderr(self, job):
        job["stderr"].extend(bytes(job["process"].readAllStandardError()))
        if len(job["stderr"]) > 64000:
            del job["stderr"][:-64000]

    def _process_error(self, job, error):
        if error == QProcess.FailedToStart:
            job["stderr"].extend(job["process"].errorString().encode("utf-8"))
            self._finished(job, -1, QProcess.CrashExit)

    def _finished(self, job, exit_code, exit_status):
        if job["generation"] not in self._jobs:
            return
        self._read_stdout(job)
        self._read_stderr(job)
        current = job["generation"] == self._generation and not job["cancelled"] and not self._closed
        if current:
            try:
                if exit_code != 0 or exit_status != QProcess.NormalExit:
                    message = job["stderr"].decode("utf-8", errors="replace").strip()
                    raise RuntimeError(message or "Geometry extraction exited unsuccessfully.")
                scene = load_scene(job["output"])
                self._publish_scene(scene)
            except Exception as error:
                self._mark_stale("OUT OF DATE — " + str(error))
        if job["cancelled"]:
            self._force_kill(job)
        self._jobs.pop(job["generation"], None)
        shutil.rmtree(job["directory"], ignore_errors=True)
        job["process"].deleteLater()
        self.cancel_button.setEnabled(any(not item["cancelled"] for item in self._jobs.values()))

    def cancel(self):
        self.debounce.stop()
        self._generation += 1
        self._cancel_jobs()
        self._mark_stale("OUT OF DATE — Update cancelled. Refresh to try again.")

    def _cancel_jobs(self):
        for job in list(self._jobs.values()):
            if not job["cancelled"]:
                job["cancelled"] = True
                self._stop_job(job)
        self.cancel_button.setEnabled(False)

    def _stop_job(self, job):
        process = job["process"]
        self._check_group(job)
        sent = False
        if job["group_owned"]:
            try:
                os.killpg(job["pid"], signal.SIGTERM)
                sent = True
            except ProcessLookupError:
                pass
        if not sent and process.state() != QProcess.NotRunning:
            process.terminate()
        QTimer.singleShot(1200, lambda: self._force_kill(job) if job["generation"] in self._jobs else None)

    def _force_kill(self, job):
        if job["group_owned"]:
            try:
                os.killpg(job["pid"], signal.SIGKILL)
            except ProcessLookupError:
                pass
        if job["process"].state() != QProcess.NotRunning:
            job["process"].kill()

    def _watch(self, files, replace_existing=False):
        if replace_existing:
            self._watched_files = set(files)
        else:
            self._watched_files.update(files)
        desired_files = {str(Path(path).resolve()) for path in self._watched_files}
        self._watched_files = desired_files
        self._file_signatures = {path: self._file_signature(path) for path in desired_files}
        desired_dirs = {str(Path(path).parent) for path in desired_files if Path(path).parent.is_dir()}
        existing_files = {path for path in desired_files if Path(path).is_file()}
        current_files = set(self.watcher.files())
        current_dirs = set(self.watcher.directories())
        remove = (current_files - existing_files) | (current_dirs - desired_dirs)
        if remove:
            self.watcher.removePaths(sorted(remove))
        add = (existing_files - current_files) | (desired_dirs - current_dirs)
        if add:
            self.watcher.addPaths(sorted(add))

    def _file_changed(self, path):
        self._file_signatures[path] = self._file_signature(path)
        self._watch([path])
        self.source_changed()

    def _directory_changed(self, path):
        missing = self._watched_files - set(self.watcher.files())
        if any(str(Path(item).parent) == path and Path(item).is_file() for item in missing):
            self._watch([])
            self.source_changed()

    @staticmethod
    def _file_signature(path):
        try:
            metadata = Path(path).stat()
            return metadata.st_mtime_ns, metadata.st_size, metadata.st_ino
        except OSError:
            return None

    def _poll_files(self):
        signatures = {path: self._file_signature(path) for path in self._watched_files}
        if signatures != self._file_signatures:
            self._file_signatures = signatures
            self._watch([])
            self.source_changed()

    def _publish_scene(self, scene):
        if any(self._file_signature(path) != signature for path, signature in self._file_signatures.items()):
            self._poll_files()
            return
        previous = self.scene
        self.scene = scene
        try:
            mesh = build_mesh(
                scene, repeats=(self.repeat_x.value(), self.repeat_y.value()),
                slice_axis=self._slice_axis_value(), slice_fraction=self.slice_slider.value() / 100,
                max_faces=30000, smooth_boundaries=self.smooth_cylinders.isChecked(),
            )
        except Exception:
            self.scene = previous
            raise
        self.mesh = mesh
        self.smooth_cylinders.setVisible(any(layer.radial is not None for layer in scene.layers))
        files = set(scene.files)
        if self._last_request is not None:
            files.update(request_files(self._last_request))
        self._watch(files, replace_existing=True)
        if not self._color_selected and any(layer.continuous for layer in scene.layers):
            self.color_mode.blockSignals(True)
            self.color_mode.setCurrentIndex(1)
            self.color_mode.blockSignals(False)
        self.repeat_y_label.setText("Invariant y extent" if scene.dimension == 1 else "Repeat a₂")
        self.repeat_y.setToolTip("Finite displayed length only; the 1D grating is infinite along y." if scene.dimension == 1 else "Repeated cells along the second lattice vector.")
        self._populate_materials()
        self.stale = False
        self.stale_artist.set_visible(False)
        self.export_button.setEnabled(True)
        if self._render():
            self.status_label.setText("Current geometry — λ = {:.8g} µm".format(scene.wavelength_um))

    def _slice_axis_value(self):
        value = self.slice_axis.currentText()
        return None if value == "None" else value.lower()

    def _slice_changed(self, *_):
        self.slice_slider.setEnabled(self.slice_axis.currentIndex() != 0)
        self.slice_label.setText("{}% retained".format(self.slice_slider.value()))
        self._rebuild_mesh()

    def _rebuild_mesh(self, *_):
        if self.scene is None:
            return
        try:
            mesh = build_mesh(
                self.scene, repeats=(self.repeat_x.value(), self.repeat_y.value()),
                slice_axis=self._slice_axis_value(), slice_fraction=self.slice_slider.value() / 100,
                max_faces=30000, smooth_boundaries=self.smooth_cylinders.isChecked(),
            )
            self.mesh = mesh
            self._populate_materials()
            self._render()
        except Exception as error:
            self._mark_stale("OUT OF DATE — " + str(error))

    def _material_key(self, label, tensor):
        return (label, tuple(np.asarray(tensor, dtype=complex).ravel()))

    @staticmethod
    def _material_color(index):
        return (0.8, 0.85, 0.9, 1.0) if index == 0 else colormaps["tab10"](((index - 1) % 10) / 9)

    def _populate_materials(self):
        self.material_table.blockSignals(True)
        try:
            for row in range(self.material_table.rowCount()):
                previous = self.material_table.cellWidget(row, 1)
                if previous is not None:
                    previous.hide()
            self.material_table.setRowCount(len(self.mesh.materials))
            self._table_keys = []
            for row, (label, tensor) in enumerate(self.mesh.materials):
                key = self._material_key(label, tensor)
                self._table_keys.append(key)
                visible, alpha = self._material_state.setdefault(key, (row != 0, 0.25 if row == 0 else 1.0))
                name_item = QTableWidgetItem(label)
                color = self._material_color(row)
                name_item.setData(Qt.DecorationRole, QColor.fromRgbF(*color))
                name_item.setFlags(Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsUserCheckable)
                name_item.setCheckState(Qt.Checked if visible else Qt.Unchecked)
                self.material_table.setItem(row, 0, name_item)
                opacity = QDoubleSpinBox()
                opacity.setRange(0, 1)
                opacity.setSingleStep(0.1)
                opacity.setDecimals(2)
                opacity.setValue(alpha)
                opacity.valueChanged.connect(lambda value, index=row: self._opacity_changed(index, value))
                self.material_table.setCellWidget(row, 1, opacity)
                tensor = np.asarray(tensor).ravel()
                details = "ε = ({})   μ = ({})".format(
                    ", ".join(self._complex_text(value) for value in tensor[:3]),
                    ", ".join(self._complex_text(value) for value in tensor[3:6]),
                )
                tensor_item = QTableWidgetItem(details)
                tensor_item.setFlags(Qt.ItemIsEnabled | Qt.ItemIsSelectable)
                tensor_item.setToolTip(details + "\nComponent order: x, y, z; tensors use native grating axes.")
                name_item.setToolTip(details)
                self.material_table.setItem(row, 2, tensor_item)
        finally:
            self.material_table.blockSignals(False)

    @staticmethod
    def _complex_text(value):
        value = complex(value)
        if value.imag == 0:
            return "{:.5g}".format(value.real)
        return "{:.5g}{:+.5g}i".format(value.real, value.imag)

    def _material_changed(self, item):
        if item.column() == 0 and item.row() < len(self._table_keys):
            key = self._table_keys[item.row()]
            self._material_state[key] = (item.checkState() == Qt.Checked, self._material_state[key][1])
            self._render()

    def _opacity_changed(self, row, value):
        if row < len(self._table_keys):
            key = self._table_keys[row]
            self._material_state[key] = (self._material_state[key][0], value)
            self._render()

    def _color_changed(self, *_):
        self._color_selected = True
        self._render()

    def _camera(self):
        return (self.axes.elev, self.axes.azim, getattr(self.axes, "roll", 0),
                self.axes.get_xlim3d(), self.axes.get_ylim3d(), self.axes.get_zlim3d())

    def _render(self):
        if self.scene is None or self.mesh is None:
            return False
        try:
            self._render_scene()
            return True
        except Exception as error:
            self._mark_stale("OUT OF DATE — Rendering failed: " + str(error))
            return False

    def _render_scene(self):
        camera = self._camera()
        if self._colorbar is not None:
            self._colorbar.remove()
            self._colorbar = None
        self.axes.clear()
        self.axes.set_position([0.10, 0.08, 0.80, 0.81])
        self.axes.set_anchor("C")
        rgba = np.empty((len(self.mesh.faces), 4), dtype=float)
        visible = np.zeros(len(self.mesh.faces), dtype=bool)
        for index, (label, tensor) in enumerate(self.mesh.materials):
            key = self._material_key(label, tensor)
            enabled, alpha = self._material_state.get(key, (index != 0, 0.9))
            match = self.mesh.material_ids == index
            rgba[match] = self._material_color(index)
            rgba[match, 3] = alpha
            visible[match] = enabled and alpha > 0
        mode = self.color_mode.currentData()
        if mode is not None and np.any(visible):
            component, part = mode
            complex_values = self.mesh.tensors[:, component]
            values = complex_values.real if part == "real" else complex_values.imag
            finite = values[visible & np.isfinite(values)]
            if finite.size:
                low, high = float(np.min(finite)), float(np.max(finite))
                if low == high:
                    margin = max(abs(low) * 0.01, 0.01)
                    low, high = low - margin, high + margin
                normalization = Normalize(low, high)
                scalar = ScalarMappable(norm=normalization, cmap="viridis")
                alpha = rgba[:, 3].copy()
                rgba[:] = scalar.to_rgba(values)
                rgba[:, 3] = alpha
                colorbar_axes = self.figure.add_axes([0.92, 0.22, 0.018, 0.55])
                self._colorbar = self.figure.colorbar(scalar, cax=colorbar_axes)
                self._colorbar.set_label(self.color_mode.currentText())
        if np.any(visible):
            faces = Poly3DCollection(
                self.mesh.faces[visible], facecolors=rgba[visible], linewidths=0, antialiased=False,
                shade=mode is None, lightsource=LightSource(azdeg=315, altdeg=45), zorder=2,
            )
            faces.set_edgecolors("none")
            self.axes.add_collection3d(faces)
        for line in self.mesh.lattice_lines:
            points = np.asarray(line)
            self.axes.plot(points[:, 0], points[:, 1], points[:, 2], color="#3e4852", linewidth=0.7, alpha=0.7, zorder=3)
        start, direction, length = self._incident_ray()
        self.incident_beam = self.axes.quiver(
            *start, *direction, length=length, color="#c63b26",
            arrow_length_ratio=0.22, linewidth=2.8, zorder=20,
        )
        self.incident_beam.set_path_effects([patheffects.Stroke(linewidth=4.4, foreground="white"), patheffects.Normal()])
        target = start + length * direction
        self.axes.scatter(*target, s=30, color="#c63b26", edgecolors="white", linewidths=0.8, depthshade=False, zorder=21)
        self.axes.text(*start, " Incident beam", color="#b42f20", fontsize=10, zorder=22,
                       bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.85, "pad": 2})
        self.axes.set_xlabel("x (µm)", labelpad=8)
        self.axes.set_ylabel("y (µm)", labelpad=8)
        self.axes.set_zlabel("z (µm)", labelpad=8)
        self.axes.set_title("{} · λ = {:.7g} µm".format(self.scene.model, self.scene.wavelength_um), pad=10)
        if self._camera_initialized:
            self.axes.view_init(elev=camera[0], azim=camera[1], roll=camera[2])
            self.axes.set_xlim3d(camera[3])
            self.axes.set_ylim3d(camera[4])
            self.axes.set_zlim3d(camera[5])
            self._set_aspect()
        else:
            self.axes.view_init(elev=26, azim=-58)
            self.fit(draw=False)
            self._camera_initialized = True
        self._update_axis_visibility()
        notes = list(dict.fromkeys(self.scene.notes + self.mesh.notes))
        lattice = np.asarray(self.scene.lattice)
        dimensions = "a₁ = {:.6g} µm".format(np.linalg.norm(lattice[0]))
        if self.scene.dimension == 2:
            dimensions += "; a₂ = {:.6g} µm".format(np.linalg.norm(lattice[1]))
        else:
            dimensions += "; displayed invariant y length = {:.6g} µm (finite illustration)".format(np.ptp(self.mesh.bounds[1]))
        dimensions += "; stack thickness = {:.6g} µm; {:,} displayed faces".format(self.scene.thickness, len(self.mesh.faces))
        if any(layer.radial is not None for layer in self.scene.layers):
            dimensions += "; smooth cylinder contours" if self.smooth_cylinders.isChecked() else "; sampled grid"
        self.resolution_label.setText(dimensions + ("\n" + " · ".join(notes) if notes else ""))
        self.canvas.draw_idle()

    def _incident_ray(self):
        bounds = np.asarray(self.mesh.bounds)
        length = max(float(np.max(np.ptp(bounds, axis=1))) * 0.25, 0.05)
        direction = incident_direction(self.scene)
        target = np.array([np.mean(bounds[0]), np.mean(bounds[1]), -self.scene.thickness if self.scene.configuration >= 2 else 0.0])
        return target - length * direction, direction, length

    def _set_aspect(self):
        extents = np.array([
            np.ptp(self.axes.get_xlim3d()), np.ptp(self.axes.get_ylim3d()), np.ptp(self.axes.get_zlim3d()),
        ])
        self.axes.set_box_aspect(np.maximum(extents, np.max(extents) * 0.05))

    def fit(self, *_args, draw=True):
        if self.mesh is None:
            return
        bounds = np.array(self.mesh.bounds, dtype=float, copy=True)
        start, _, _ = self._incident_ray()
        bounds[:, 0] = np.minimum(bounds[:, 0], start)
        bounds[:, 1] = np.maximum(bounds[:, 1], start)
        extent = np.maximum(np.ptp(bounds, axis=1), 0.01)
        bounds[:, 0] -= extent * 0.08
        bounds[:, 1] += extent * 0.08
        self.axes.set_xlim3d(bounds[0])
        self.axes.set_ylim3d(bounds[1])
        self.axes.set_zlim3d(bounds[2])
        self._set_aspect()
        if draw:
            self.canvas.draw_idle()

    def reset_view(self):
        self.axes.view_init(elev=26, azim=-58)
        self._update_axis_visibility()
        self.fit()

    def set_view(self, elevation, azimuth):
        self.axes.view_init(elev=elevation, azim=azimuth)
        self._update_axis_visibility()
        self.canvas.draw_idle()

    def _update_axis_visibility(self):
        elevation = math.radians(self.axes.elev)
        azimuth = math.radians(self.axes.azim)
        eye = np.array([math.cos(elevation) * math.cos(azimuth), math.cos(elevation) * math.sin(azimuth), math.sin(elevation)])
        for axis, alignment in zip((self.axes.xaxis, self.axes.yaxis, self.axes.zaxis), np.abs(eye)):
            visible = alignment < .999999
            axis.set_visible(visible)
            axis.label.set_visible(visible)
            axis.offsetText.set_visible(visible)
            axis.line.set_visible(visible)
            axis.gridlines.set_visible(visible)
            axis.set_major_locator(AutoLocator() if visible else NullLocator())

    def _camera_changed(self, event):
        if event.inaxes == self.axes:
            self._update_axis_visibility()
            self.canvas.draw_idle()

    def _scroll(self, event):
        if event.inaxes != self.axes:
            return
        scale = 0.85 if event.button == "up" else 1 / 0.85
        for getter, setter in ((self.axes.get_xlim3d, self.axes.set_xlim3d), (self.axes.get_ylim3d, self.axes.set_ylim3d), (self.axes.get_zlim3d, self.axes.set_zlim3d)):
            low, high = getter()
            middle = (low + high) / 2
            half = (high - low) * scale / 2
            setter(middle - half, middle + half)
        self.canvas.draw_idle()

    def export_png(self, path=None):
        if self.scene is None:
            return False
        if path is None:
            path, _ = QFileDialog.getSaveFileName(self, "Export 3D Geometry", "rcw-geometry.png", "PNG image (*.png)")
        if not path:
            return False
        export_artists = []
        try:
            if self.color_mode.currentData() is None:
                handles = []
                for index, (label, tensor) in enumerate(self.mesh.materials):
                    visible, opacity = self._material_state[self._material_key(label, tensor)]
                    if visible and opacity:
                        handles.append(Patch(facecolor=self._material_color(index), alpha=opacity, label=label))
                if handles:
                    export_artists.append(self.figure.legend(
                        handles=handles[:12], loc="upper left", bbox_to_anchor=(0.02, 0.91), fontsize=8,
                        title="Materials" if len(handles) <= 12 else "Materials (first 12)",
                        framealpha=0.9,
                    ))
            export_artists.append(self.figure.text(
                0.01, 0.015, self.resolution_label.text().split("\n")[0], ha="left", va="bottom", fontsize=7,
            ))
            self.figure.savefig(path, dpi=180, facecolor="white")
            self.status_label.setText(("OUT OF DATE — " if self.stale else "") + "Exported " + str(path))
            return True
        except Exception as error:
            self.status_label.setText(("OUT OF DATE — " if self.stale else "") + "Export failed: " + str(error))
            return False
        finally:
            for artist in export_artists:
                artist.remove()
            self.canvas.draw_idle()
