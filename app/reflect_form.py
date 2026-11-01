import os
import csv
import datetime
import importlib, importlib.util, sys
import shutil
import math
from pathlib import Path
from gui_validation import finite_float, optical_constant
from gui_solver import SolverCancelled, run_solver
from scatmech_paths import get_data_dir, get_scatmech_bin
from PyQt5.QtWidgets import (
    QWidget, QLabel, QLineEdit,
    QVBoxLayout, QHBoxLayout,
    QPushButton, QComboBox,
    QTextEdit, QGroupBox,
    QFormLayout,
    QTableWidget, QTableWidgetItem,
    QSizePolicy
)
from PyQt5.QtCore import QUrl
from PyQt5.QtGui import QDesktopServices
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

_scatmech_bin_path = get_scatmech_bin()
scatmech_bin = str(_scatmech_bin_path)
if _scatmech_bin_path.is_dir():
    os.environ["PATH"] = scatmech_bin + os.pathsep + os.environ.get("PATH", "")


def parse_reflect_output(output):
    rows = []
    for line_number, line in enumerate(output.splitlines(), start=1):
        parts = line.split()
        if not parts:
            continue
        if not rows and [part.lower() for part in parts] == ["theta", "rp", "rs"]:
            continue
        if len(parts) != 3:
            raise ValueError(f"Reflectance output row {line_number} has {len(parts)} columns; expected 3.")
        try:
            values = [float(part) for part in parts]
        except ValueError as exc:
            raise ValueError(f"Reflectance output row {line_number} contains a nonnumeric value.") from exc
        if not all(math.isfinite(value) for value in values):
            raise ValueError(f"Reflectance output row {line_number} contains a nonfinite value.")
        rows.append(parts)
    if not rows:
        raise ValueError("No finite reflectance rows found in stdout.")
    return rows


class ReflectForm(QWidget):

    def __init__(self):
        super().__init__()

        self.main_layout = QHBoxLayout(self)
        self.setLayout(self.main_layout)


        self.form_layout = QVBoxLayout()
        form_widget = QWidget()
        form_widget.setLayout(self.form_layout)


        self.figure = Figure(figsize=(5.6, 3.4), tight_layout=True)
        self.canvas = FigureCanvas(self.figure)
        self.canvas.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)

        self.main_layout.addWidget(form_widget, 0)
        self.main_layout.addWidget(self.canvas, 1)


        wl_group = QGroupBox("Wavelength [µm]")
        wl_form = QFormLayout()
        self.wavelength = QLineEdit("0.633")
        wl_form.addRow("λ (µm):", self.wavelength)
        wl_group.setLayout(wl_form)


        sub_group = QGroupBox("Substrate dielectric (n,k)")
        sub_form = QFormLayout()
        self.sub_n = QLineEdit("1.50")
        self.sub_k = QLineEdit("0.00")
        row_subnk = QHBoxLayout()
        row_subnk.addWidget(QLabel("n:"))
        row_subnk.addWidget(self.sub_n)
        row_subnk.addWidget(QLabel("k:"))
        row_subnk.addWidget(self.sub_k)
        sub_form.addRow("Substrate (n,k):", row_subnk)
        sub_group.setLayout(sub_form)


        stack_group = QGroupBox("Film Stack (bottom → top)")
        stack_v = QVBoxLayout()

        tool_row = QHBoxLayout()
        self.btn_add = QPushButton("Add Layer")
        self.btn_remove = QPushButton("Remove Selected")
        self.btn_up = QPushButton("Move Up")
        self.btn_down = QPushButton("Move Down")
        for b in (self.btn_add, self.btn_remove, self.btn_up, self.btn_down):
            tool_row.addWidget(b)
        tool_row.addStretch(1)
        stack_v.addLayout(tool_row)

        self.tbl = QTableWidget(0, 2)
        self.tbl.setHorizontalHeaderLabels(["Material (n,k) or file", "Thickness [µm]"])
        self.tbl.setToolTip("Layers are listed from substrate to surface. Use (n,k), a real index, or a dielectric file.")
        self.tbl.horizontalHeader().setStretchLastSection(True)
        stack_v.addWidget(self.tbl)

        stack_group.setLayout(stack_v)


        plot_group = QGroupBox("Plot")
        plot_form = QFormLayout()
        self.plot_column = QComboBox()
        self.plot_column.addItems(["R_p (p-pol)", "R_s (s-pol)"])
        plot_form.addRow("Y-axis:", self.plot_column)
        plot_group.setLayout(plot_form)


        action_row = QHBoxLayout()
        self.run_btn = QPushButton("Run reflectprog")
        self.clear_btn = QPushButton("Clear Plot")
        self.open_output_btn = QPushButton("Open Last Output")
        self.open_input_btn = QPushButton("Open Last Input")
        action_row.addWidget(self.run_btn, 1)
        action_row.addWidget(self.clear_btn, 1)
        action_row.addWidget(self.open_input_btn, 1)
        action_row.addWidget(self.open_output_btn, 1)


        self.output_box = QTextEdit()
        self.output_box.setReadOnly(True)
        self.output_box.setMinimumHeight(120)


        self.form_layout.addWidget(wl_group)
        self.form_layout.addWidget(sub_group)
        self.form_layout.addWidget(stack_group)
        self.form_layout.addWidget(plot_group)
        self.form_layout.addLayout(action_row)
        self.form_layout.addWidget(QLabel("Log:"))
        self.form_layout.addWidget(self.output_box)
        self.form_layout.addStretch(1)


        self.btn_add.clicked.connect(lambda checked=False: self._add_layer())
        self.btn_remove.clicked.connect(self._remove_selected)
        self.btn_up.clicked.connect(lambda: self._move_selected(-1))
        self.btn_down.clicked.connect(lambda: self._move_selected(+1))

        self.run_btn.clicked.connect(self.run_reflectprog)
        self.clear_btn.clicked.connect(self.clear_plot)
        self.open_output_btn.clicked.connect(self.open_last_output)
        self.open_input_btn.clicked.connect(self.open_last_input)


        self.last_stdout_path = None
        self.last_input_path = None
        self.last_csv_path = None
        self.last_plot_success = False
        self._running = False
        self.plot_column.currentIndexChanged.connect(self._update_plot_selection)


        self._add_layer(default_material="(1.50,0.00)", default_thickness="0.100000")



    def _add_layer(self, default_material=None, default_thickness=None):
        r = self.tbl.rowCount()
        self.tbl.insertRow(r)
        mat = QTableWidgetItem(str(default_material) if default_material is not None else "(1.50,0.00)")
        thk = QTableWidgetItem(str(default_thickness) if default_thickness is not None else "0.100000")
        self.tbl.setItem(r, 0, mat)
        self.tbl.setItem(r, 1, thk)

    def _remove_selected(self):
        rows = {idx.row() for idx in self.tbl.selectedIndexes()}
        for r in sorted(rows, reverse=True):
            self.tbl.removeRow(r)

    def _move_selected(self, direction: int):
        r = self.tbl.currentRow()
        if r < 0:
            return
        new_r = r + direction
        if not (0 <= new_r < self.tbl.rowCount()):
            return
        row_vals = [self.tbl.item(r, c).text() if self.tbl.item(r, c) else "" for c in range(2)]
        self.tbl.removeRow(r)
        self.tbl.insertRow(new_r)
        self.tbl.setItem(new_r, 0, QTableWidgetItem(row_vals[0]))
        self.tbl.setItem(new_r, 1, QTableWidgetItem(row_vals[1]))
        self.tbl.selectRow(new_r)


    def clear_plot(self):
        self.last_plot_success = False
        self.last_csv_path = None
        self.figure.clear()
        self.canvas.draw()

    def open_last_output(self):
        self._open_file(self.last_stdout_path, "output")

    def open_last_input(self):
        self._open_file(self.last_input_path, "input")

    def _open_file(self, path, description):
        try:
            if not path or not Path(path).is_file():
                self.output_box.append(f"No {description} file to open.")
                return
            if not QDesktopServices.openUrl(QUrl.fromLocalFile(str(Path(path).resolve()))):
                self.output_box.append(f"Could not open {description} file: {path}")
        except OSError as exc:
            self.output_box.append(f"Could not open {description} file: {exc}")

    def _update_plot_selection(self):
        if self.last_csv_path:
            self.render_with_external(self.last_csv_path, self.plot_column.currentIndex() + 1)

    def to_params(self):
        return {
            "wavelength_um": self.wavelength.text().strip(),
            "substrate": {"n": self.sub_n.text().strip(), "k": self.sub_k.text().strip()},
            "layers": [
                {
                    "material": self.tbl.item(row, 0).text().strip() if self.tbl.item(row, 0) else "",
                    "thickness_um": self.tbl.item(row, 1).text().strip() if self.tbl.item(row, 1) else "",
                }
                for row in range(self.tbl.rowCount())
            ],
            "plot_y": self.plot_column.currentText(),
        }

    def from_params(self, params):
        if not isinstance(params, dict):
            raise ValueError("Reflect parameters must be a dictionary.")
        candidate = self.to_params()
        wavelength = params.get("wavelength_um", params.get("wavelength"))
        if wavelength is not None:
            candidate["wavelength_um"] = wavelength
        substrate = params.get("substrate")
        if substrate is None:
            substrate = {key: params[f"sub_{key}"] for key in ("n", "k") if f"sub_{key}" in params}
        if isinstance(substrate, dict):
            candidate["substrate"].update(substrate)
        else:
            normalized = optical_constant(substrate, "Substrate")
            if not normalized.startswith("("):
                raise ValueError("Reflect substrate must be a real index or an (n,k) pair.")
            n, k = normalized[1:-1].split(",")
            candidate["substrate"] = {"n": n, "k": k}

        if "layers" in params:
            if not isinstance(params["layers"], list):
                raise ValueError("Reflect layers must be a list.")
            layers = []
            for index, layer in enumerate(params["layers"], start=1):
                if not isinstance(layer, dict):
                    raise ValueError(f"Layer {index} must be a dictionary.")
                material = layer.get("material")
                if material is None:
                    n = finite_float(layer.get("n", 1.5), f"Layer {index} n")
                    k = finite_float(layer.get("k", 0), f"Layer {index} k")
                    material = f"({n},{k})"
                layers.append({
                    "material": material,
                    "thickness_um": layer.get("thickness_um", layer.get("thickness", "0.1")),
                })
            candidate["layers"] = layers

        plot_y = params.get("plot_y", candidate["plot_y"])
        aliases = {"p": "R_p (p-pol)", "rp": "R_p (p-pol)", "r_p": "R_p (p-pol)",
                   "s": "R_s (s-pol)", "rs": "R_s (s-pol)", "r_s": "R_s (s-pol)"}
        if isinstance(plot_y, str):
            plot_y = aliases.get(plot_y.strip().lower(), plot_y.strip())
        if plot_y not in ("R_p (p-pol)", "R_s (s-pol)"):
            raise ValueError("Reflect plot_y must select R_p (p-pol) or R_s (s-pol).")
        candidate["plot_y"] = plot_y
        self._parameter_input_lines(candidate)

        self.wavelength.setText(str(candidate["wavelength_um"]))
        self.sub_n.setText(str(candidate["substrate"]["n"]))
        self.sub_k.setText(str(candidate["substrate"]["k"]))
        self.tbl.setRowCount(0)
        for layer in candidate["layers"]:
            self._add_layer(layer["material"], layer["thickness_um"])
        previous = self.plot_column.blockSignals(True)
        self.plot_column.setCurrentText(plot_y)
        self.plot_column.blockSignals(previous)


    def _build_input_lines(self):
        return self._parameter_input_lines(self.to_params())

    @staticmethod
    def _parameter_input_lines(params):
        wavelength = finite_float(params["wavelength_um"], "Wavelength (µm)", minimum=0, strict_minimum=True)
        n = finite_float(params["substrate"]["n"], "Substrate n")
        k = finite_float(params["substrate"]["k"], "Substrate k")
        substrate = optical_constant(f"({n},{k})", "Substrate")


        layer_pairs = []
        for r, layer in enumerate(params["layers"]):
            mat_txt = layer["material"]
            thk_txt = layer["thickness_um"]
            if mat_txt == "" and thk_txt == "":
                continue
            material = optical_constant(mat_txt, f"Layer {r+1} material")
            thickness = finite_float(thk_txt, f"Layer {r+1} thickness (µm)", minimum=0)
            layer_pairs.extend([material, str(thickness)])

        lines = []
        lines.append(str(wavelength))
        lines.append(substrate)
        lines.append(" ".join(layer_pairs) if layer_pairs else "")
        lines.append("")
        return lines


    def render_with_external(self, csv_path: str, y_idx: int):
        self.last_plot_success = False
        self.figure.clear()
        ax = self.figure.gca()

        here = Path(__file__).resolve().parent
        cwd = Path(os.getcwd())
        csv_dir = Path(csv_path).resolve().parent if csv_path else None

        tried = []
        mod = None
        how = None

        def _try_import(name):
            nonlocal mod, how
            try:
                mod = importlib.import_module(name)
                how = f"import {name}"
                return True
            except Exception as e:
                tried.append(f"import {name}: {e}")
                return False

        def _try_path(path):
            nonlocal mod, how
            try:
                if path and path.exists():
                    spec = importlib.util.spec_from_file_location("reflectplot", str(path))
                    if spec and spec.loader:
                        mod = importlib.util.module_from_spec(spec)
                        sys.modules["reflectplot"] = mod
                        spec.loader.exec_module(mod)
                        how = f"spec_from_file_location({path})"
                        return True
            except Exception as e:
                tried.append(f"load {path}: {e}")
            return False

        if not _try_import("reflectplot"):
            if not _try_path(here / "reflectplot.py"):
                if not _try_path(cwd / "reflectplot.py") and csv_dir:
                    _try_path(csv_dir / "reflectplot.py")

        if mod is None:
            self.output_box.append("Could not import reflectplot.py: " + " | ".join(tried))
            self.canvas.draw()
            return
        else:
            self.output_box.append(f"reflectplot resolved via: {how}")

        fn = getattr(mod, "plot_csv", None)
        if not callable(fn):
            self.output_box.append("reflectplot.py found, but it must define plot_csv(ax, csv_path, ...).")
            self.canvas.draw()
            return

        try:
            label = "Rp" if y_idx == 1 else ("Rs" if y_idx == 2 else None)
            fn(ax, csv_path, x_col=0, y_col=y_idx, semilogy=False, label=label)
            self.canvas.draw()
            self.last_plot_success = True
            self.output_box.append("Plot updated via reflectplot.plot_csv")
        except Exception as e:
            self.output_box.append(f"reflectplot render error: {e}")
            self.canvas.draw()


    def run_reflectprog(self):
        if self._running:
            return
        self.output_box.clear()
        self.clear_plot()
        self.last_input_path = None
        self.last_stdout_path = None
        self.last_csv_path = None
        try:
            input_lines = self._build_input_lines()
        except ValueError as exc:
            self.output_box.append(f"[Error] Invalid input: {exc}")
            return

        exe = shutil.which("reflectprog")
        if not exe:
            self.output_box.setText(
                "[Error] 'reflectprog' not found. Ensure it is on PATH or next to the app.\n"
                f"Current PATH includes: {scatmech_bin}"
            )
            return

        timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S_%f")
        input_text = "\n".join(input_lines) + "\n"
        try:
            data_dir = get_data_dir()
            data_dir.mkdir(parents=True, exist_ok=True)
            input_filename = data_dir / f"reflect_input_{timestamp}.txt"
            output_filename = data_dir / f"reflect_output_{timestamp}.txt"
            csv_filename = data_dir / f"reflect_output_{timestamp}.csv"
            input_filename.write_text(input_text, encoding="utf-8")
            self.last_input_path = str(input_filename)
        except OSError as exc:
            self.output_box.append(f"[Error] Could not save reflectprog input: {exc}")
            return

        self._running = True
        self.run_btn.setEnabled(False)
        try:
            result = run_solver(exe, input=input_text, parent=self, cwd=data_dir)
        except SolverCancelled:
            self.output_box.append("reflectprog cancelled.")
            return
        except (OSError, RuntimeError) as exc:
            self.output_box.append(f"[Error] Could not invoke reflectprog: {exc}")
            return
        finally:
            self._running = False
            self.run_btn.setEnabled(True)

        stdout = result.stdout or ""
        stderr = result.stderr or ""
        try:
            output_filename.write_text(
                stdout + ("\n--- STDERR ---\n" + stderr if stderr else ""), encoding="utf-8"
            )
            self.last_stdout_path = str(output_filename)
            self.output_box.append(f"Saved output log: {output_filename}")
        except OSError as exc:
            self.output_box.append(f"[Error] Could not save reflectprog output: {exc}")
            return

        if result.returncode != 0:
            self.output_box.append(f"[Error] reflectprog failed (exit code {result.returncode}):\n{stderr or stdout}")
            return

        if "SCATMECH: Error" in stderr:
            self.output_box.append(f"[Error] reflectprog reported a SCATMECH error:\n{stderr}")
            return
        try:
            rows = parse_reflect_output(stdout)
        except ValueError as exc:
            self.output_box.append(f"[Error] {exc}")
            if stderr:
                self.output_box.append(stderr)
            return

        try:
            with csv_filename.open("w", newline="", encoding="utf-8") as csvfile:
                writer = csv.writer(csvfile)
                writer.writerow(["theta", "Rp", "Rs"])
                writer.writerows(rows)
        except OSError as exc:
            self.output_box.append(f"[Error] Could not save reflectance CSV: {exc}")
            return

        self.last_csv_path = str(csv_filename)
        self.output_box.append(f"Saved CSV: {csv_filename}")
        self._update_plot_selection()
