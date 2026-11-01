import os
import csv
import datetime
import math
import shutil
from pathlib import Path
from scatmech_paths import get_scatmech_bin, get_data_dir
from gui_solver import SolverCancelled, run_solver
from gui_validation import finite_float, optical_constant
from PyQt5.QtWidgets import (
    QWidget, QLabel, QLineEdit,
    QVBoxLayout, QHBoxLayout,
    QPushButton, QTextEdit, QGroupBox,
    QFormLayout,
    QDialog, QTableWidget, QTableWidgetItem, QSizePolicy, QMenu
)
from PyQt5.QtCore import pyqtSignal

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

_scatmech_bin_path = get_scatmech_bin()
scatmech_bin = str(_scatmech_bin_path)
if _scatmech_bin_path.is_dir():
    os.environ["PATH"] = scatmech_bin + os.pathsep + os.environ.get("PATH", "")


def parse_mie_output(output):

    lines = output.splitlines()
    for index, line in enumerate(lines):
        header = line.split()
        if header and header[0].lower() in {"theta", "angle"}:
            break
    else:
        raise ValueError("Could not find a Mie data header ('Theta' or 'Angle').")
    names = [name.lower() for name in header]
    if len(set(names)) != len(names):
        raise ValueError("Mie output contains duplicate column names.")
    if not any(name in names for name in ("s11", "pol", "s33", "s34")):
        raise ValueError("Mie output has no recognized scattering columns.")
    rows = []
    for line in lines[index + 1:]:
        parts = line.split()
        if not parts:
            continue
        if len(parts) != len(header):
            raise ValueError(f"Mie output row {len(rows) + 1} has {len(parts)} columns; expected {len(header)}.")
        try:
            values = [float(value) for value in parts]
        except ValueError as exc:
            raise ValueError(f"Mie output row {len(rows) + 1} contains a nonnumeric value.") from exc
        for column, name in enumerate(names):
            if name in {"theta", "angle", "phi"} and not math.isfinite(values[column]):
                raise ValueError("Mie output contains a nonfinite scattering angle.")
        rows.append(parts)
    if not rows:
        raise ValueError("No data rows found in Mie output.")
    return header, rows


class MieForm(QWidget):

    requestClearPlot = pyqtSignal()

    def __init__(self):
        super().__init__()


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


        grp = QGroupBox("Parameters")
        form = QFormLayout()

        self.step_angle_deg = QLineEdit("1")
        form.addRow("Step Angle (theta, deg):", self.step_angle_deg)

        self.step_phi_deg = QLineEdit("1")
        form.addRow("Step Azimuth (phi, deg):", self.step_phi_deg)

        self.wavelength_um = QLineEdit("0.532")
        form.addRow("Wavelength (µm):", self.wavelength_um)

        self.medium_optics = QLineEdit("(1,0)")
        form.addRow("Nonabsorbing Surrounding Medium (n,0):", self.medium_optics)
        self.medium_optics.setToolTip("Use a positive real index, an (n,0) pair, or a nonabsorbing dielectric file.")

        self.radius_um = QLineEdit("0.05")
        form.addRow("Radius (µm):", self.radius_um)

        self.sphere_optics = QLineEdit("(1.59,0)")
        form.addRow("Optical Properties of the Sphere (n,k):", self.sphere_optics)

        grp.setLayout(form)
        self.form_layout.addWidget(grp)


        ctrl = QHBoxLayout()
        self.run_btn = QPushButton("Run MieProg")
        self.metric_menu = QMenu(self.run_btn)
        for _name in ["S11", "Pol", "S33", "S34"]:
            act = self.metric_menu.addAction(_name)
            act.triggered.connect(lambda checked, n=_name: self.run_with_metric(n))

        self.run_btn.setMenu(self.metric_menu)

        self.run_btn.setStyleSheet(
    "QPushButton { text-align: center; } "
    "QPushButton::menu-indicator { subcontrol-position: right center; }"
)

        self.clear_btn = QPushButton("Clear Plot")
        self.log_color_btn = QPushButton("Log Color: Off")
        self.log_color_btn.setCheckable(True)
        self.open_output_btn = QPushButton("Open Last Output")
        self.open_input_btn = QPushButton("Open Last Input")

        ctrl.setContentsMargins(0, 0, 0, 0)
        ctrl.setSpacing(12)

        ctrl.addWidget(self.run_btn,        1)
        ctrl.addWidget(self.clear_btn,      1)
        ctrl.addWidget(self.log_color_btn,  1)
        ctrl.addWidget(self.open_output_btn,1)
        ctrl.addWidget(self.open_input_btn, 1)

        ctrl_wrap = QWidget()
        ctrl_wrap.setLayout(ctrl)
        self.form_layout.addWidget(ctrl_wrap)
        self.form_layout.addWidget(QLabel("Log:"))


        self.log = QTextEdit()
        self.log.setReadOnly(True)
        self.log.setMinimumHeight(140)
        self.log.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.form_layout.addWidget(self.log)
        self.form_layout.addStretch(1)


        self.run_btn.clicked.connect(self.run_mieprog)
        self.clear_btn.clicked.connect(self.clear_plot)
        self.log_color_btn.toggled.connect(self.toggle_log_color_scale)
        self.open_output_btn.clicked.connect(self.open_last_output)
        self.open_input_btn.clicked.connect(self.open_last_input)


        self._plot_clear_callback = None


        self.metric_name = "S11"
        self.color_scale = "linear"


        self.last_stdout_path = None
        self.last_input_path = None
        self.last_csv_path = None
        self.last_plot_success = False
        self._running = False

    def connect_plot_clear(self, slot):
        self._plot_clear_callback = slot


    def clear_plot(self):
        self.last_plot_success = False
        self.last_csv_path = None
        sig = getattr(self, "requestClearPlot", None)
        if sig is not None:
            try:
                sig.emit()
            except Exception as e:
                self.log.append(f"Clear-plot signal error: {e}")
        if callable(getattr(self, "_plot_clear_callback", None)):
            try:
                self._plot_clear_callback()
            except Exception as e:
                self.log.append(f"Clear-plot callback error: {e}")
        if hasattr(self, 'figure'):
            self.figure.clear()
        if hasattr(self, 'canvas'):
            self.canvas.draw()
        self.log.append("Plot cleared.")

    def run_mieprog(self):
        if self._running:
            return
        self._running = True
        self.run_btn.setEnabled(False)
        self.last_plot_success = False
        self.last_input_path = None
        self.last_stdout_path = None
        self.last_csv_path = None
        try:
            self.figure.clear()
            self.canvas.draw()
            stdin_payload = self._input_payload()
            exe = shutil.which("mieprog")
            if not exe:
                self.log.append(
                    "[Error] 'mieprog' not found. Set SCATMECH_BIN or add its directory to PATH.\n"
                    f"Configured solver directory: {get_scatmech_bin()}"
                )
                return
            timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S_%f")
            data_dir = get_data_dir()
            data_dir.mkdir(parents=True, exist_ok=True)
            input_txt = str(data_dir / f"mie_input_{timestamp}.txt")
            stdout_txt = str(data_dir / f"mie_output_{timestamp}.txt")
            csv_filename = str(data_dir / f"mie_output_{timestamp}.csv")
            Path(input_txt).write_text(stdin_payload, encoding="utf-8")
            self.last_input_path = input_txt
            self.log.append(f"Saved input deck: {input_txt}")
            self.log.append(f"Running MieProg: {exe}")
            proc = run_solver(exe, input=stdin_payload, parent=self,
                              cwd=data_dir)

            with open(stdout_txt, "w", encoding="utf-8", errors="ignore") as f:
                f.write(proc.stdout or "")
                if proc.stderr:
                    f.write("\n--- STDERR ---\n")
                    f.write(proc.stderr)
            self.last_stdout_path = stdout_txt

            if proc.returncode != 0:
                self.log.append(f"mieprog exited with code {proc.returncode}. See output text for details.")
                return
            if "SCATMECH: Error" in (proc.stderr or ""):
                self.log.append("mieprog reported a SCATMECH error. See output text for details.")
                return

            self.log.append("mieprog completed. Parsing output table…")
            header, data_rows = parse_mie_output(proc.stdout or "")
            names = [name.upper() for name in header]
            if self.metric_name.upper() not in names:
                raise ValueError(f"Mie output does not contain the selected {self.metric_name} metric.")
            metric_column = names.index(self.metric_name.upper())
            if not any(math.isfinite(float(row[metric_column])) for row in data_rows):
                raise ValueError(f"Mie output has no finite {self.metric_name} values.")


            with open(csv_filename, "w", newline="", encoding="utf-8") as fcsv:
                w = csv.writer(fcsv)
                w.writerow(header)
                w.writerows(data_rows)
            self.last_csv_path = csv_filename
            self.log.append(f"Saved CSV: {csv_filename}")

            self.render_with_external(csv_filename)

        except SolverCancelled:
            self.log.append("mieprog cancelled.")
        except ValueError as e:
            self.log.append(f"Mie input/output error: {e}")
        except OSError as e:
            self.log.append(f"Mie file or executable error: {e}")
        except Exception as e:
            self.log.append(f"Error running mieprog: {e}")
        finally:
            self._running = False
            self.run_btn.setEnabled(True)

    def _input_payload(self):
        params = self.to_params()
        numeric = {}
        for field, label in (
            ("step_angle_deg", "Theta step"), ("step_phi_deg", "Phi step"),
            ("wavelength_um", "Wavelength"), ("radius_um", "Sphere radius"),
        ):
            numeric[field] = finite_float(params[field], label, minimum=0, strict_minimum=True)

        for field, label, endpoint in (
            ("step_angle_deg", "Theta", 180), ("step_phi_deg", "Phi", 360),
        ):
            bound = math.radians(endpoint)
            if bound + math.radians(numeric[field]) <= bound:
                raise ValueError(f"{label} step is too small to advance the angular scan.")
        if not isinstance(self.metric_name, str) or self.metric_name not in {"S11", "Pol", "S33", "S34"}:
            raise ValueError("Mie metric must be S11, Pol, S33, or S34.")
        if self.color_scale not in {"linear", "log"}:
            raise ValueError("Mie color scale must be linear or log.")
        medium = optical_constant(params["medium_optics"], "Surrounding medium")
        if medium.startswith("("):
            n, k = (float(value) for value in medium[1:-1].split(","))
            if n <= 0 or k != 0:
                raise ValueError("Mie surrounding medium must have a positive index and zero absorption (k=0).")
        values = [
            str(numeric["step_angle_deg"]), str(numeric["step_phi_deg"]),
            str(numeric["wavelength_um"]),
            medium,
            str(numeric["radius_um"]),
            optical_constant(params["sphere_optics"], "Sphere optical properties"),
        ]
        return "\n".join(values) + "\n"


    def render_with_external(self, csv_path: str):
        self.last_plot_success = False
        try:
            import mieplot

            mieplot.set_metric(self.metric_name)
            mieplot.set_color_scale(self.color_scale)
            self.figure.clear()
            ax = self.figure.add_subplot(111)
            mieplot.plot_csv(ax, csv_path)
            self.canvas.draw()
            self.last_plot_success = True
            self.log.append("Plot updated via mieplot.plot_csv")
            return True
        except Exception as exc:
            self.log.append(f"mieplot render error: {exc}")
            self.canvas.draw()
            return False

    def toggle_log_color_scale(self, checked: bool):
        self.color_scale = "log" if checked else "linear"
        self.log_color_btn.setText("Log Color: On" if checked else "Log Color: Off")
        csv_path = getattr(self, "last_csv_path", None)
        if not csv_path or not os.path.exists(csv_path):
            view_name = "log" if checked else "normal"
            self.log.append(f"Color scale selected: {view_name}. It will apply after the next MieProg run.")
            return

        view_name = "log" if checked else "normal"
        self.log.append(f"Color scale selected: {view_name}. Re-rendering last Mie plot...")
        self.render_with_external(csv_path)

    def run_with_metric(self, name: str):
        if self._running:
            return
        if not isinstance(name, str) or name not in {"S11", "Pol", "S33", "S34"}:
            self.log.append("Mie metric must be S11, Pol, S33, or S34.")
            return
        self.metric_name = name
        self.log.append(f"Metric selected: {name}. Running mieprog...")
        self.run_mieprog()



    def to_params(self):
        return {
            "step_angle_deg": self._get_text("step_angle_deg", "1"),
            "step_phi_deg": self._get_text("step_phi_deg", "10"),
            "wavelength_um": self._get_text("wavelength_um", "0.532"),
            "medium_optics": self._get_text("medium_optics", "(1,0)"),
            "radius_um": self._get_text("radius_um", "0.05"),
            "sphere_optics": self._get_text("sphere_optics", "(1.59,0)"),
            "metric": self.metric_name,
            "color_scale": self.color_scale,
        }

    def from_params(self, p: dict):
        if not isinstance(p, dict):
            raise ValueError("Mie parameters must be a dictionary.")
        metric = p.get("metric", p.get("metric_name", self.metric_name))
        scale = str(p.get("color_scale", self.color_scale)).strip().lower()
        scale = {"normal": "linear", "lin": "linear", "logarithmic": "log"}.get(scale, scale)
        if not isinstance(metric, str) or metric not in {"S11", "Pol", "S33", "S34"}:
            raise ValueError("Mie metric must be S11, Pol, S33, or S34.")
        if scale not in {"linear", "log"}:
            raise ValueError("Mie color scale must be linear or log.")
        self._set_text("step_angle_deg", p.get("step_angle_deg"))
        self._set_text("step_phi_deg", p.get("step_phi_deg"))
        self._set_text("wavelength_um", p.get("wavelength_um"))
        self._set_text("medium_optics", p.get("medium_optics"))
        self._set_text("radius_um", p.get("radius_um"))
        self._set_text("sphere_optics", p.get("sphere_optics"))
        self.metric_name = metric
        self.color_scale = scale
        previous = self.log_color_btn.blockSignals(True)
        self.log_color_btn.setChecked(scale == "log")
        self.log_color_btn.setText("Log Color: On" if scale == "log" else "Log Color: Off")
        self.log_color_btn.blockSignals(previous)


    def _find_latest(self, prefix: str, folder: str):
        try:
            paths = [os.path.join(folder, f) for f in os.listdir(folder)
                     if f.startswith(prefix) and Path(folder, f).is_file()
                     and Path(f).suffix.lower() in {".txt", ".csv"}]
            if not paths:
                return None
            return max(paths, key=lambda p: os.path.getmtime(p))
        except Exception:
            return None

    def _read_file(self, path: str):
        try:
            with open(path, "r", encoding="utf-8", errors="ignore") as f:
                return f.read()
        except Exception as e:
            return f"(Could not open file: {e})"

    def open_last_output(self):
        import csv as _csv
        data_dir = get_data_dir()
        path = getattr(self, "last_stdout_path", None)
        if not path or not os.path.exists(path):
            path = self._find_latest("mie_output_", data_dir)
        if not path:
            self.log.append(f"No Mie output file found in {data_dir}.")
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

        _, ext = os.path.splitext(path)
        ext = (ext or "").lower()
        if ext == ".csv":
            try:
                with open(path, newline="") as f:
                    rows = list(_csv.reader(f))
                if not rows:
                    txt.setPlainText("(Empty CSV)")
                else:
                    header, data = rows[0], rows[1:]
                    if any(len(row) != len(header) for row in data):
                        raise ValueError("CSV data rows do not match the header column count.")
                    table.setColumnCount(len(header))
                    table.setRowCount(len(data))
                    table.setHorizontalHeaderLabels(header)
                    for r, row in enumerate(data):
                        for c, val in enumerate(row):
                            table.setItem(r, c, QTableWidgetItem(val))
                    table.setSortingEnabled(True)
                    table.resizeColumnsToContents()
                    txt.setPlainText("")
            except Exception as e:
                txt.setPlainText(f"(Failed to open CSV: {e})")
        else:
            text = self._read_file(path)
            stdout, _, stderr = text.partition("\n--- STDERR ---\n")
            lines = stdout.splitlines()
            header_idx = None
            for i, line in enumerate(lines):
                parts = line.split()
                if parts and parts[0].lower() in {"theta", "angle"}:
                    header_idx = i
                    break
            pre = "\n".join(lines[:header_idx]) if header_idx not in (None, 0) else ""
            if stderr:
                pre += "\n--- STDERR ---\n" + stderr
            txt.setPlainText(pre if pre else "(No preamble)")
            if header_idx is not None:
                try:
                    header, data_rows = parse_mie_output(stdout)
                except ValueError as exc:
                    txt.setPlainText(f"{text}\n\nCould not display table: {exc}")
                else:
                    table.setColumnCount(len(header))
                    table.setRowCount(len(data_rows))
                    table.setHorizontalHeaderLabels(header)
                    for r, row in enumerate(data_rows):
                        for c, val in enumerate(row):
                            table.setItem(r, c, QTableWidgetItem(val))
                    table.setSortingEnabled(True)
                    table.resizeColumnsToContents()
            else:
                txt.setPlainText(text)

        dlg.resize(900, 700)
        dlg.exec_()

    def open_last_input(self):
        data_dir = get_data_dir()
        path = getattr(self, "last_input_path", None)
        if not path or not os.path.exists(path):
            path = self._find_latest("mie_input_", data_dir)
        if not path:
            self.log.append(f"No Mie input deck found in {data_dir}.")
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


    def _get_text(self, name, default=""):
        try:
            w = getattr(self, name, None)
            return w.text().strip() if w is not None else default
        except Exception:
            return default

    def _set_text(self, name, value):
        try:
            w = getattr(self, name, None)
            if w is not None and value is not None:
                w.setText(str(value))
        except Exception:
            pass
