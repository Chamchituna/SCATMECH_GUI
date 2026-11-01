import os
import re
import csv
import datetime
import math
import importlib
import importlib.util
import sys
import shutil
from pathlib import Path
from scatmech_paths import get_scatmech_bin, get_data_dir
from scatmech_gratings import (
    ONE_D_GRATING_SPECS, get_one_d_grating_spec,
    serialize_one_d_grating, validate_one_d_grating,
)
from gui_validation import finite_float, integer, single_line
from gui_solver import SolverCancelled, run_solver
from custom_profile_binding import CustomProfileBinding, prepare_profiles
from PyQt5.QtWidgets import (
    QWidget, QLabel, QLineEdit,
    QVBoxLayout, QHBoxLayout,
    QPushButton, QTextEdit, QGroupBox,
    QFormLayout, QFileDialog, QComboBox,
    QSizePolicy, QTableWidget, QTableWidgetItem, QDialog
)
from PyQt5.QtCore import pyqtSignal

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

_scatmech_bin_path = get_scatmech_bin()
scatmech_bin = str(_scatmech_bin_path)
if _scatmech_bin_path.is_dir():
    os.environ["PATH"] = scatmech_bin + os.pathsep + os.environ.get("PATH", "")


_TYPE_CHOICES = [
    ("Forward / Reflection", "0"),
    ("Forward / Transmission", "1"),
    ("Backward / Reflection", "2"),
    ("Backward / Transmission", "3"),
]

class RCWForm(QWidget):

    requestClearPlot = pyqtSignal()

    def __init__(self):
        super().__init__()

        self._suspend_preview_log = False
        self._geometry_viewer = None
        self._grating_param_store = {}
        self._current_grating_module = None
        self.grating_param_inputs = {}

        self.main_layout = QHBoxLayout(self)
        self.setLayout(self.main_layout)

        self.form_layout = QVBoxLayout()
        form_widget = QWidget()
        form_widget.setLayout(self.form_layout)

        self.figure = Figure(figsize=(6, 5), tight_layout=True)
        self.canvas = FigureCanvas(self.figure)
        self.canvas.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)

        plot_layout = QVBoxLayout()
        plot_layout.addWidget(self.canvas, 1)
        plot_widget = QWidget()
        plot_widget.setLayout(plot_layout)

        self.main_layout.addWidget(form_widget, 1)
        self.main_layout.addWidget(plot_widget, 1)
        self.main_layout.setStretch(0, 1)
        self.main_layout.setStretch(1, 1)

        basic_group = QGroupBox("Simulation Parameters")
        basic_form = QFormLayout()

        self.order = QLineEdit("6")
        basic_form.addRow("Maximum Diffraction Order:", self.order)

        self.type_combo = QComboBox()
        for label, code in _TYPE_CHOICES:
            self.type_combo.addItem(label, code)
        self.type_combo.setToolTip(
            "Transmission and backward configurations require a nonabsorbing transmission medium (k=0)."
        )
        basic_form.addRow("Configuration Type:", self.type_combo)

        self.wavelength_um = QLineEdit("0.532")
        basic_form.addRow("Initial Model Wavelength (µm):", self.wavelength_um)
        self.wavelength_um.setToolTip("Initial model value. Each calculation uses the wavelength sweep below.")

        self.theta_inc_deg = QLineEdit("0")
        basic_form.addRow("Incident Polar Angle θᵢ (deg):", self.theta_inc_deg)

        self.rotation_deg = QLineEdit("0")
        basic_form.addRow("Grating Rotation φ (deg):", self.rotation_deg)

        basic_group.setLayout(basic_form)
        self.form_layout.addWidget(basic_group)

        sweep_group = QGroupBox("Wavelength Sweep")
        sweep_form = QFormLayout(sweep_group)
        self.sweep_start_um = QLineEdit("0.2")
        self.sweep_end_um = QLineEdit("0.8")
        self.sweep_step_um = QLineEdit("0.005")
        sweep_form.addRow("Start (µm):", self.sweep_start_um)
        sweep_form.addRow("End (µm):", self.sweep_end_um)
        sweep_form.addRow("Step (µm):", self.sweep_step_um)
        self.form_layout.addWidget(sweep_group)

        grating_group = QGroupBox("Grating Definition")
        grating_form = QFormLayout()

        self.grating_model_combo = QComboBox()
        for model_name in ONE_D_GRATING_SPECS:
            self.grating_model_combo.addItem(model_name)
        grating_form.addRow("Module:", self.grating_model_combo)

        self.grating_param_container = QWidget()
        self.grating_param_layout = QFormLayout()
        self.grating_param_container.setLayout(self.grating_param_layout)
        grating_form.addRow(self.grating_param_container)

        self.browse_grating = None
        self.custom_profile_binding = CustomProfileBinding(self, 1, self._custom_profile_node,
            self._set_custom_profile_node, self._custom_profile_context, self._lock_custom_profile_fields,
            self._on_grating_param_changed)
        self.custom_profile_btn = self.custom_profile_binding.button
        grating_form.addRow(self.custom_profile_binding.controls)
        grating_group.setLayout(grating_form)
        self.form_layout.addWidget(grating_group)

        ctrl_row = QHBoxLayout()
        self.run_btn = QPushButton("Run RCWProg")
        self.clear_btn = QPushButton("Clear Plot")
        self.open_output_btn = QPushButton("Open Last Output")
        self.open_input_btn = QPushButton("Open Last Input")
        ctrl_row.addWidget(self.run_btn, 1)
        ctrl_row.addWidget(self.clear_btn, 1)
        ctrl_row.addWidget(self.open_input_btn, 1)
        ctrl_row.addWidget(self.open_output_btn, 1)
        ctrl_widget = QWidget()
        ctrl_widget.setLayout(ctrl_row)
        self.form_layout.addWidget(ctrl_widget)
        self.geometry_btn = QPushButton("View 3D Geometry")
        self.geometry_btn.clicked.connect(self.show_geometry_viewer)
        self.form_layout.addWidget(self.geometry_btn)

        self.form_layout.addWidget(QLabel("Input Deck Preview:"))
        self.input_preview = QTextEdit()
        self.input_preview.setReadOnly(True)
        self.input_preview.setMinimumHeight(140)
        self.form_layout.addWidget(self.input_preview)

        self.form_layout.addWidget(QLabel("Log:"))
        self.log = QTextEdit()
        self.log.setReadOnly(True)
        self.log.setMinimumHeight(140)
        self.log.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.form_layout.addWidget(self.log)
        self.form_layout.addStretch(1)

        self.run_btn.clicked.connect(self.run_rcwprog)
        self.clear_btn.clicked.connect(self.clear_plot)
        self.open_output_btn.clicked.connect(self.open_last_output)
        self.open_input_btn.clicked.connect(self.open_last_input)
        self.grating_model_combo.currentTextChanged.connect(self._on_grating_model_changed)
        for field in (self.order, self.wavelength_um, self.theta_inc_deg, self.rotation_deg,
                      self.sweep_start_um, self.sweep_end_um, self.sweep_step_um):
            field.textChanged.connect(self._on_grating_param_changed)
        self.type_combo.currentIndexChanged.connect(self._on_grating_param_changed)

        self._plot_clear_callback = None
        self.last_stdout_path = None
        self.last_input_path = None
        self.last_csv_path = None
        self.last_plot_success = False

        self._rebuild_grating_param_form(self.grating_model_combo.currentText())
        self.populate_input_preview()

    def connect_plot_clear(self, slot):
        self._plot_clear_callback = slot

    def _custom_profile_node(self):
        model = self.grating_model_combo.currentText()
        return {"kind": "grating", "model": model, "params": {
            field["name"]: self._get_grating_param_value(field["name"], field["default"])
            for field in ONE_D_GRATING_SPECS[model]["fields"]}}

    def _set_custom_profile_node(self, node):
        self._grating_param_store[node["model"]] = dict(node["params"])
        blocked = self.grating_model_combo.blockSignals(True)
        try:
            self.grating_model_combo.setCurrentText(node["model"])
        finally:
            self.grating_model_combo.blockSignals(blocked)
        self._rebuild_grating_param_form(node["model"])

    def _custom_profile_context(self):
        return {"wavelength_um": finite_float(self.sweep_start_um.text(), "Sweep start wavelength", minimum=0, strict_minimum=True),
                "incident_deg": finite_float(self.theta_inc_deg.text(), "Incident polar angle", minimum=-90, maximum=90),
                "rotation_deg": finite_float(self.rotation_deg.text(), "Grating rotation"),
                "configuration": int(self.type_combo.currentData())}

    def _lock_custom_profile_fields(self, active):
        for widget in self.grating_param_inputs.values():
            widget.setReadOnly(active)
        if self.browse_grating is not None:
            self.browse_grating.setEnabled(not active)

    def geometry_request(self):
        self.custom_profile_binding.ensure()
        from geometry_types import GeometryRequest

        model = self.grating_model_combo.currentText()
        if model not in ONE_D_GRATING_SPECS:
            raise ValueError("Select a supported RCW grating.")
        node = {"kind": "grating", "model": model, "params": {
            field["name"]: self._get_grating_param_value(field["name"], field["default"])
            for field in ONE_D_GRATING_SPECS[model]["fields"]
        }}
        errors = validate_one_d_grating(node)
        if errors:
            raise ValueError("\n".join(errors))
        wavelength = finite_float(self.sweep_start_um.text(), "Sweep start wavelength", minimum=0, strict_minimum=True)
        incident = finite_float(self.theta_inc_deg.text(), "Incident polar angle", minimum=-90, maximum=90)
        if abs(incident) == 90:
            raise ValueError("Incident polar angle must be strictly between -90 and 90 degrees.")
        return GeometryRequest(
            grating=node, wavelength_um=wavelength, incident_deg=incident,
            rotation_deg=finite_float(self.rotation_deg.text(), "Grating rotation"),
            configuration=int(self.type_combo.currentData()), source_label="RCW",
        )

    def show_geometry_viewer(self):
        if self._geometry_viewer is None:
            from geometry_viewer import GeometryViewer
            self._geometry_viewer = GeometryViewer(parent=self, request_provider=self.geometry_request)
        self._geometry_viewer.show()
        self._geometry_viewer.raise_()
        self._geometry_viewer.activateWindow()

    def _geometry_source_changed(self, *_):
        if self._geometry_viewer is not None:
            self._geometry_viewer.source_changed()

    def closeEvent(self, event):
        self.custom_profile_binding.close()
        if self._geometry_viewer is not None:
            self._geometry_viewer.close()
        super().closeEvent(event)

    def clear_plot(self):
        self.last_plot_success = False
        sig = getattr(self, "requestClearPlot", None)
        if sig is not None:
            try:
                sig.emit()
            except Exception as exc:
                self.log.append(f"Clear-plot signal error: {exc}")
        if callable(getattr(self, "_plot_clear_callback", None)):
            try:
                self._plot_clear_callback()
            except Exception as exc:
                self.log.append(f"Clear-plot callback error: {exc}")
        if hasattr(self, "figure"):
            self.figure.clear()
        if hasattr(self, "canvas"):
            self.canvas.draw()
        self.log.append("Plot cleared.")

    def _browse_for_grating(self):
        fname, _ = QFileDialog.getOpenFileName(self, "Select grating file", "", "All Files (*)")
        if fname:
            widget = self.grating_param_inputs.get("filename")
            if widget is not None:
                widget.setText(fname)
            else:
                self.log.append("Filename field not available for current grating module.")

    def populate_input_preview(self):
        try:
            payload = self._build_input_payload()
        except ValueError as exc:
            payload = f"Invalid input: {exc}"
        self.input_preview.setPlainText(payload)
        if not getattr(self, "_suspend_preview_log", False):
            self.log.append("Parameters updated.")

    def _build_input_payload(self):
        self.custom_profile_binding.ensure()
        order = str(integer(self.order.text(), "Maximum diffraction order", minimum=0))
        type_code = self.type_combo.currentData() or "0"
        wavelength = self.wavelength_um.text().strip()
        finite_float(wavelength, "Wavelength", minimum=0, strict_minimum=True)
        theta = self.theta_inc_deg.text().strip()
        angle = finite_float(theta, "Incident polar angle", minimum=-90, maximum=90)
        if abs(angle) == 90:
            raise ValueError("Incident polar angle must be strictly between -90 and 90 degrees.")
        rotation = self.rotation_deg.text().strip()
        finite_float(rotation, "Grating rotation")
        module_name = self.grating_model_combo.currentText() or next(iter(ONE_D_GRATING_SPECS))
        node = {"model": module_name, "params": {
            field["name"]: self._get_grating_param_value(field["name"], field["default"])
            for field in self._get_model_fields(module_name)
        }}
        errors = validate_one_d_grating(node)
        if errors:
            raise ValueError("\n".join(errors))
        start = finite_float(self.sweep_start_um.text(), "Sweep start", minimum=0, strict_minimum=True)
        end = finite_float(self.sweep_end_um.text(), "Sweep end", minimum=start)
        step = finite_float(self.sweep_step_um.text(), "Sweep step", minimum=0, strict_minimum=True)
        if start + step == start or end + step == end:
            raise ValueError("Sweep step is too small to advance the wavelength.")

        lines = [order, type_code, wavelength, theta, rotation]
        lines.extend(serialize_one_d_grating(node))
        lines.extend(["lambda", self.sweep_start_um.text().strip(),
                      self.sweep_end_um.text().strip(), self.sweep_step_um.text().strip()])

        return "\n".join(lines) + "\n"

    def _get_model_fields(self, module_name):
        if module_name not in ONE_D_GRATING_SPECS:
            module_name = next(iter(ONE_D_GRATING_SPECS))
        return get_one_d_grating_spec(module_name)["fields"]

    def _get_grating_param_value(self, name, default):
        widget = self.grating_param_inputs.get(name)
        if widget is None:
            return default
        return widget.text().strip()

    def _store_current_grating_params(self):
        module = getattr(self, "_current_grating_module", None)
        if not module:
            return
        values = {}
        for name, widget in self.grating_param_inputs.items():
            values[name] = widget.text()
        self._grating_param_store[module] = values

    def _rebuild_grating_param_form(self, module_name):
        while self.grating_param_layout.count():
            item = self.grating_param_layout.takeAt(0)
            widget = item.widget()
            if widget is not None:
                widget.hide()
                widget.deleteLater()
        self.grating_param_inputs = {}
        self.browse_grating = None

        values = self._grating_param_store.get(module_name, {})

        for field in self._get_model_fields(module_name):
            name = field["name"]
            stored_value = values.get(name)
            value = stored_value if stored_value is not None else field["default"]
            use_browse = module_name == "Generic_Grating" and name == "filename"
            self._add_grating_param_row(name, field["label"], value, browse=use_browse)

        self._current_grating_module = module_name

    def _add_grating_param_row(self, name, label, value, browse=False):
        line_edit = QLineEdit()
        line_edit.setText(value or "")
        line_edit.textChanged.connect(self._on_grating_param_changed)
        self.grating_param_inputs[name] = line_edit

        if browse:
            row_widget = QWidget()
            row_layout = QHBoxLayout()
            row_layout.setContentsMargins(0, 0, 0, 0)
            row_layout.addWidget(line_edit)
            browse_btn = QPushButton("Browse…")
            browse_btn.clicked.connect(self._browse_for_grating)
            row_layout.addWidget(browse_btn)
            row_widget.setLayout(row_layout)
            self.grating_param_layout.addRow(f"{label}:", row_widget)
            self.browse_grating = browse_btn
        else:
            self.grating_param_layout.addRow(f"{label}:", line_edit)

    def _on_grating_model_changed(self, module_name):
        self.custom_profile_binding.detach()
        self._store_current_grating_params()
        self._rebuild_grating_param_form(module_name)
        self._on_grating_param_changed(module_name)

    def _on_grating_param_changed(self, _value=None):
        previous_state = self._suspend_preview_log
        self._suspend_preview_log = True
        try:
            self.populate_input_preview()
        finally:
            self._suspend_preview_log = previous_state
        self._geometry_source_changed()

    def run_rcwprog(self):
        if getattr(self, "_solver_running", False):
            self.log.append("A calculation is already running for this module.")
            return
        self.last_plot_success = False
        self.last_csv_path = None
        self.last_stdout_path = None
        self.last_input_path = None
        self.figure.clear()
        self.canvas.draw()
        try:
            payload = self._build_input_payload()
        except ValueError as exc:
            self.log.append(f"[Input error] {exc}")
            return
        self.input_preview.setPlainText(payload)

        timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S_%f")
        data_dir = get_data_dir()
        input_txt = os.path.join(data_dir, f"rcw_input_{timestamp}.txt")
        output_txt = os.path.join(data_dir, f"rcw_stdout_{timestamp}.log")
        csv_filename = os.path.join(data_dir, f"rcw_output_{timestamp}.csv")

        try:
            os.makedirs(data_dir, exist_ok=True)
            with open(input_txt, "w", encoding="utf-8") as f:
                f.write(payload)
            self.last_input_path = input_txt
            self.log.append(f"Saved input deck: {input_txt}")
        except Exception as exc:
            self.log.append(f"[Error] Could not write input deck: {exc}")
            return

        exe = shutil.which("rcwprog")
        if not exe:
            self.log.append(
                "[Error] 'rcwprog' not found. Ensure it is on PATH or next to the app.\n"
                f"Current PATH includes: {scatmech_bin}"
            )
            return

        self.log.append(f"Running RCWProg: {exe}")
        try:
            proc = run_solver(
                exe,
                input=payload,
                parent=self,
                cwd=data_dir,
            )
        except SolverCancelled:
            self.log.append("rcwprog cancelled.")
            return
        except Exception as exc:
            self.log.append(f"[Error] Could not invoke rcwprog: {exc}")
            return

        try:
            with open(output_txt, "w", encoding="utf-8", errors="ignore") as f:
                f.write(proc.stdout or "")
                if proc.stderr:
                    f.write("\n--- STDERR ---\n")
                    f.write(proc.stderr)
            self.last_stdout_path = output_txt
            self.log.append(f"Saved stdout log: {output_txt}")
        except Exception as exc:
            self.log.append(f"[Warning] Failed to store stdout: {exc}")

        if proc.returncode != 0 or "SCATMECH: Error" in (proc.stderr or "") + (proc.stdout or ""):
            self.log.append("rcwprog failed. See output log for details.")
            self.log.append(proc.stderr or proc.stdout or "")
            return

        self.log.append("rcwprog completed. Parsing output table…")

        try:
            header, rows = self._extract_table(proc.stdout)
        except Exception as exc:
            self.log.append(f"[Error] Failed to parse rcwprog output: {exc}")
            return

        if not rows:
            self.log.append("No numeric rows detected in rcwprog output.")
            self.log.append(proc.stderr or "")
            return

        try:
            with open(csv_filename, "w", newline="", encoding="utf-8") as fcsv:
                writer = csv.writer(fcsv)
                if header:
                    writer.writerow(header)
                writer.writerows(rows)
            self.log.append(f"Saved CSV: {csv_filename}")
            self.last_csv_path = csv_filename
        except Exception as exc:
            self.log.append(f"[Error] Could not write CSV: {exc}")
            return

        self.render_with_external(csv_filename)

    def _extract_table(self, stdout: str):
        if not stdout:
            return [], []
        lines = stdout.splitlines()
        header_idx = None
        header = []
        for i, line in enumerate(lines):
            tokens = re.split(r"\s+", line.strip())
            if not tokens:
                continue
            lowered = [t.lower() for t in tokens]
            if len(lowered) >= 2 and (
                ({"alpha", "beta"}.issubset(lowered) or {"α", "β"}.issubset(lowered))
                or lowered[0] in {"lambda", "wavelength", "λ", "order", "theta", "phi"}
            ):
                header_idx = i
                header = tokens
                break
        data_rows = []
        start = header_idx + 1 if header_idx is not None else 0
        for line in lines[start:]:
            s = line.strip()
            if not s:
                continue
            parts = re.split(r"\s+", s)
            if not parts:
                continue
            numeric_prefix = 0
            for token in parts:
                try:
                    float(token)
                    numeric_prefix += 1
                except Exception:
                    break
            if numeric_prefix == 0:
                if data_rows:
                    raise ValueError("RCW output contains a nonnumeric table row.")
                continue
            if numeric_prefix != len(parts) or (header and len(parts) != len(header)):
                raise ValueError("RCW output contains an incomplete or nonnumeric table row.")
            if not math.isfinite(float(parts[0])):
                raise ValueError("RCW output contains a nonfinite scan coordinate.")
            data_rows.append(parts)
        if not data_rows:
            return header, []
        if not header:
            width = max(len(row) for row in data_rows)
            header = [f"col{i+1}" for i in range(width)]
        if any(len(row) != len(header) for row in data_rows):
            raise ValueError("RCW output contains inconsistent column counts.")
        return header, data_rows

    def render_with_external(self, csv_path: str):
        self.last_plot_success = False
        self.figure.clear()
        ax = self.figure.gca()

        here = os.path.dirname(os.path.abspath(__file__))
        cwd = os.getcwd()
        csv_dir = os.path.dirname(os.path.abspath(csv_path)) if csv_path else None

        tried = []

        def _try_normal(name):
            try:
                mod = importlib.import_module(name)
                return mod, f"import {name}"
            except Exception as exc:
                tried.append(f"import {name}: {exc}")
                return None, None

        def _try_file(path):
            try:
                if path and os.path.exists(path):
                    spec = importlib.util.spec_from_file_location("rcwplot", path)
                    if spec and spec.loader:
                        mod = importlib.util.module_from_spec(spec)
                        sys.modules["rcwplot"] = mod
                        spec.loader.exec_module(mod)
                        return mod, f"load {path}"
            except Exception as exc:
                tried.append(f"load {path}: {exc}")
            return None, None

        if here not in sys.path:
            sys.path.insert(0, here)

        mod, how = _try_normal("rcwplot")
        if mod is None:
            mod, how = _try_file(os.path.join(here, "rcwplot.py"))
        if mod is None:
            mod, how = _try_file(os.path.join(cwd, "rcwplot.py"))
        if mod is None and csv_dir:
            mod, how = _try_file(os.path.join(csv_dir, "rcwplot.py"))

        if mod is None:
            self.log.append("Could not import rcwplot.py: " + " | ".join(tried))
            self.canvas.draw()
            return
        else:
            self.log.append(f"rcwplot resolved via: {how}")

        fn = getattr(mod, "plot_csv", None)
        if not callable(fn):
            self.log.append("rcwplot.py found, but it must define plot_csv(ax, csv_path).")
            self.canvas.draw()
            return

        try:
            fn(ax, csv_path)
            self.canvas.draw()
            self.last_plot_success = True
            self.log.append("Plot updated via rcwplot.plot_csv")
        except Exception as exc:
            self.log.append(f"rcwplot render error: {exc}")
            self.canvas.draw()

    def open_last_output(self):
        data_dir = get_data_dir()
        path = getattr(self, "last_csv_path", None)
        if not path or not os.path.exists(path):
            path = self._find_latest("rcw_output_", data_dir, suffix=".csv")
        if not path:
            self.log.append(f"No RCW CSV output found in {data_dir}.")
            return

        dlg = QDialog(self)
        dlg.setWindowTitle(f"Output: {os.path.basename(path)}")
        layout = QVBoxLayout(dlg)

        txt = QTextEdit(dlg)
        txt.setReadOnly(True)
        txt.setMinimumHeight(120)
        layout.addWidget(txt)

        table = QTableWidget(dlg)
        table.setMinimumHeight(320)
        layout.addWidget(table)

        btn_row = QHBoxLayout()
        close_btn = QPushButton("Close", dlg)
        btn_row.addStretch(1)
        btn_row.addWidget(close_btn)
        layout.addLayout(btn_row)
        close_btn.clicked.connect(dlg.accept)

        try:
            with open(path, newline="") as f:
                rows = list(csv.reader(f))
            if not rows:
                txt.setPlainText("(Empty CSV)")
            else:
                header, data = rows[0], rows[1:]
                table.setColumnCount(len(header))
                table.setRowCount(len(data))
                table.setHorizontalHeaderLabels(header)
                for r, row in enumerate(data):
                    for c, val in enumerate(row):
                        table.setItem(r, c, QTableWidgetItem(val))
                table.setSortingEnabled(True)
                table.resizeColumnsToContents()
                txt_path = getattr(self, "last_stdout_path", None)
                if txt_path and os.path.exists(txt_path):
                    txt.setPlainText(self._read_file(txt_path))

                else:
                    txt.setPlainText("(Stdout log not available)")
        except Exception as exc:
            txt.setPlainText(f"(Failed to open CSV: {exc})")

        dlg.resize(900, 700)
        dlg.exec_()

    def open_last_input(self):
        data_dir = get_data_dir()
        path = getattr(self, "last_input_path", None)
        if not path or not os.path.exists(path):
            path = self._find_latest("rcw_input_", data_dir)
        if not path:
            self.log.append(f"No RCW input deck found in {data_dir}.")
            return
        content = self._read_file(path)
        dlg = QDialog(self)
        dlg.setWindowTitle(f"Input: {os.path.basename(path)}")
        layout = QVBoxLayout(dlg)
        view = QTextEdit(dlg)
        view.setReadOnly(True)
        view.setPlainText(content)
        layout.addWidget(view)
        btn_row = QHBoxLayout()
        close_btn = QPushButton("Close", dlg)
        btn_row.addStretch(1)
        btn_row.addWidget(close_btn)
        layout.addLayout(btn_row)
        close_btn.clicked.connect(dlg.accept)
        dlg.resize(700, 500)
        dlg.exec_()

    def _find_latest(self, prefix: str, folder: str, suffix: str = ""):
        try:
            paths = [
                os.path.join(folder, f)
                for f in os.listdir(folder)
                if f.startswith(prefix) and f.endswith(suffix)
            ] if suffix else [
                os.path.join(folder, f)
                for f in os.listdir(folder)
                if f.startswith(prefix)
            ]
            if not paths:
                return None
            return max(paths, key=lambda p: os.path.getmtime(p))
        except Exception:
            return None

    def _read_file(self, path: str):
        try:
            with open(path, "r", encoding="utf-8", errors="ignore") as f:
                return f.read()
        except Exception as exc:
            return f"(Could not open file: {exc})"

    def to_params(self):
        self.custom_profile_binding.ensure()
        return {
            "custom_profiles": {"grating": self.custom_profile_binding.saved()} if self.custom_profile_binding.document is not None else {},
            "order": self.order.text().strip(),
            "type": self.type_combo.currentData(),
            "wavelength_um": self.wavelength_um.text().strip(),
            "theta_inc_deg": self.theta_inc_deg.text().strip(),
            "rotation_deg": self.rotation_deg.text().strip(),
            "sweep_start_um": self.sweep_start_um.text().strip(),
            "sweep_end_um": self.sweep_end_um.text().strip(),
            "sweep_step_um": self.sweep_step_um.text().strip(),
            "grating_module": self.grating_model_combo.currentText(),
            "grating_params": {
                name: widget.text().strip()
                for name, widget in self.grating_param_inputs.items()
            },
            "input_deck": self.input_preview.toPlainText(),
        }

    def from_params(self, params: dict):
        if not isinstance(params, dict):
            raise ValueError("RCW parameters must be an object.")
        if not params:
            return
        params = dict(params)
        profile_documents = {}
        if params.get("custom_profiles"):
            if not isinstance(params.get("grating_params", {}), dict):
                raise ValueError("Grating parameters must be an object.")
            node = {"kind": "grating", "model": params.get("grating_module", self.grating_model_combo.currentText()),
                    "params": dict(params.get("grating_params", {}))}
            profile_documents, node = prepare_profiles(params["custom_profiles"], node)
            params["grating_module"] = node["model"]
            params["grating_params"] = node["params"]
        elif not isinstance(params.get("custom_profiles", {}), dict):
            raise ValueError("custom_profiles must be an object.")
        field_names = ("order", "wavelength_um", "theta_inc_deg", "rotation_deg",
                       "sweep_start_um", "sweep_end_um", "sweep_step_um")
        values = {name: single_line(params[name], name, allow_empty=True)
                  for name in field_names if name in params}
        type_code = str(params.get("type", self.type_combo.currentData()))
        type_index = self.type_combo.findData(type_code)
        if type_index < 0:
            raise ValueError("RCW configuration type must be 0, 1, 2, or 3.")
        module_name = params.get("grating_module", self.grating_model_combo.currentText())
        if not isinstance(module_name, str) or module_name not in ONE_D_GRATING_SPECS:
            raise ValueError(f"Unsupported grating module: {module_name}")
        stored_params = params.get("grating_params", {})
        if not isinstance(stored_params, dict):
            raise ValueError("Grating parameters must be an object.")
        valid_names = {field["name"] for field in self._get_model_fields(module_name)}
        unknown = set(stored_params) - valid_names
        if unknown:
            raise ValueError(f"Unknown grating parameters: {', '.join(sorted(unknown))}")
        stored_params = {name: single_line(value, name, allow_empty=True)
                         for name, value in stored_params.items()}

        previous_state = self._suspend_preview_log
        self._suspend_preview_log = True
        try:
            self.custom_profile_binding.restore(None, None)
            for name, value in values.items():
                getattr(self, name).setText(value)
            self.type_combo.setCurrentIndex(type_index)
            self._store_current_grating_params()
            updated = dict(self._grating_param_store.get(module_name, {}))
            updated.update(stored_params)
            self._grating_param_store[module_name] = updated
            blocked = self.grating_model_combo.blockSignals(True)
            try:
                self.grating_model_combo.setCurrentText(module_name)
            finally:
                self.grating_model_combo.blockSignals(blocked)
            self._rebuild_grating_param_form(module_name)
            self.custom_profile_binding.restore(profile_documents.get("grating"), self._custom_profile_node())
            self.populate_input_preview()
        finally:
            self._suspend_preview_log = previous_state
        self._geometry_source_changed()
