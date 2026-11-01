import os
import csv
import datetime
import math
import importlib, importlib.util, sys
import inspect
import shutil
from pathlib import Path
import re
import json
from scatmech_paths import get_scatmech_bin, get_data_dir
from gui_validation import finite_float, integer, optical_constant, single_line, solver_file
from gui_solver import run_solver, SolverCancelled
from custom_profile_binding import (CustomProfileBinding, editor_profiles, grating_editor_context, grating_editor_media,
    has_editor_profiles, node_expression, normalize_profile_boundaries, prepare_profiles, restore_editor_profiles)
from scatmech_gratings import (
    CROSS_GRATING_SPECS,
    ONE_D_GRATING_SPECS,
    build_default_cross_grating,
    coerce_cross_grating,
    coerce_one_d_grating,
    list_cross_grating_models,
    serialize_cross_grating,
    serialize_one_d_grating,
    validate_cross_grating,
    validate_one_d_grating,
)

from PyQt5.QtWidgets import (
    QWidget, QLabel, QLineEdit,
    QVBoxLayout, QHBoxLayout,
    QPushButton, QComboBox,
    QTextEdit, QGroupBox,
    QFormLayout, QFileDialog,
    QSizePolicy
)
from PyQt5.QtCore import pyqtSignal
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

_scatmech_bin_path = get_scatmech_bin()
scatmech_bin = str(_scatmech_bin_path)
if _scatmech_bin_path.is_dir():
    os.environ["PATH"] = scatmech_bin + os.pathsep + os.environ.get("PATH", "")


def parse_brdf_output(stdout):
    rows = []
    for line_number, raw in enumerate(stdout.splitlines(), 1):
        text = raw.strip()
        if not text or text.startswith("#"):
            continue
        if "SCATMECH: Error" in text:
            raise ValueError(text)
        tokens = text.split()
        try:
            float(tokens[0].replace("D", "E").replace("d", "e"))
        except ValueError:
            for token in tokens[1:]:
                try:
                    float(token.replace("D", "E").replace("d", "e"))
                except ValueError:
                    continue
                raise ValueError(f"BRDF result line {line_number} contains a malformed value.") from None
            continue
        if len(tokens) != 7:
            raise ValueError(f"BRDF result line {line_number} must contain seven columns.")
        try:
            row = [float(token.replace("D", "E").replace("d", "e")) for token in tokens]
        except ValueError:
            raise ValueError(f"BRDF result line {line_number} contains a malformed value.") from None
        if not all(math.isfinite(row[i]) for i in (0, 1, 6)):
            raise ValueError("brdfprog returned non-finite angles or intensity.")
        rows.append(row)
    if not rows:
        raise ValueError("brdfprog returned no result rows.")
    return rows


class OneDGratingEditor(QWidget):
    changed = pyqtSignal()

    def __init__(self, *, title: str = "1D Grating", parent=None):
        super().__init__(parent)
        self._value_store = {}
        self.param_inputs = {}
        self._current_model = ""

        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)

        self.group = QGroupBox(title)
        self.group_layout = QFormLayout()
        self.group.setLayout(self.group_layout)
        outer.addWidget(self.group)

        self.model_combo = QComboBox()
        for name in ONE_D_GRATING_SPECS:
            self.model_combo.addItem(name)
        self.group_layout.addRow("Model:", self.model_combo)

        self.fields_container = QWidget()
        self.fields_layout = QFormLayout()
        self.fields_container.setLayout(self.fields_layout)
        self.group_layout.addRow(self.fields_container)

        self.model_combo.currentTextChanged.connect(self._on_model_changed)
        initial_model = self.model_combo.currentText() or next(iter(ONE_D_GRATING_SPECS))
        self._rebuild_fields(initial_model)
        self.custom_profile_binding = CustomProfileBinding(self, 1, self.to_node, self.from_node,
            lambda: grating_editor_context(self), self._lock_custom_profile_fields, self.changed.emit,
            lambda: grating_editor_media(self))
        self.custom_profile_btn = self.custom_profile_binding.button
        self.group_layout.addRow(self.custom_profile_binding.controls)

    def _lock_custom_profile_fields(self, active):
        for widget in self.param_inputs.values():
            widget.setReadOnly(active)
        for button in self.fields_container.findChildren(QPushButton):
            if button.text() == "Browse...":
                button.setEnabled(not active)

    def closeEvent(self, event):
        self.custom_profile_binding.close()
        super().closeEvent(event)

    def _clear_fields(self):
        while self.fields_layout.count():
            item = self.fields_layout.takeAt(0)
            widget = item.widget()
            if widget is not None:
                widget.deleteLater()
        self.param_inputs = {}

    def _add_param_row(self, field, value: str):
        name = field["name"]
        line_edit = QLineEdit()
        line_edit.setText(value)
        line_edit.textChanged.connect(lambda _value: self.changed.emit())
        self.param_inputs[name] = line_edit
        if name == "filename":
            row_widget = QWidget()
            row_layout = QHBoxLayout(row_widget)
            row_layout.setContentsMargins(0, 0, 0, 0)
            row_layout.addWidget(line_edit)
            browse_btn = QPushButton("Browse...")
            browse_btn.clicked.connect(lambda: self._browse_for_file(name))
            row_layout.addWidget(browse_btn)
            self.fields_layout.addRow(f"{field['label']}:", row_widget)
        else:
            self.fields_layout.addRow(f"{field['label']}:", line_edit)

    def _browse_for_file(self, field_name: str):
        fname, _ = QFileDialog.getOpenFileName(self, "Select grating file", "", "All Files (*)")
        if fname:
            widget = self.param_inputs.get(field_name)
            if widget is not None:
                widget.setText(fname)

    def _store_current_state(self):
        if self._current_model:
            self._value_store[self._current_model] = self.to_node()

    def _rebuild_fields(self, model_name: str, node=None):
        if model_name not in ONE_D_GRATING_SPECS:
            model_name = next(iter(ONE_D_GRATING_SPECS))
        self._clear_fields()
        state = node if node is not None else self._value_store.get(model_name)
        coerced = coerce_one_d_grating(state, default_model=model_name)
        self._current_model = coerced["model"]
        self._value_store[self._current_model] = coerced
        for field in ONE_D_GRATING_SPECS[self._current_model]["fields"]:
            self._add_param_row(field, coerced["params"][field["name"]])

    def _on_model_changed(self, model_name: str):
        self.custom_profile_binding.detach()
        self._store_current_state()
        self._rebuild_fields(model_name)
        self.changed.emit()

    def to_node(self, ensure=True):
        if ensure and hasattr(self, "custom_profile_binding"):
            self.custom_profile_binding.ensure()
        model = self._current_model or self.model_combo.currentText() or next(iter(ONE_D_GRATING_SPECS))
        params = {}
        for field in ONE_D_GRATING_SPECS[model]["fields"]:
            widget = self.param_inputs.get(field["name"])
            value = widget.text().strip() if widget is not None else field["default"]
            params[field["name"]] = value
        return {
            "kind": "grating",
            "model": model,
            "params": params,
        }

    def from_node(self, node):
        coerced = coerce_one_d_grating(node)
        self._value_store[coerced["model"]] = coerced
        idx = self.model_combo.findText(coerced["model"])
        if idx >= 0:
            previous = self.model_combo.blockSignals(True)
            self.model_combo.setCurrentIndex(idx)
            self.model_combo.blockSignals(previous)
        self._rebuild_fields(coerced["model"], node=coerced)
        self.changed.emit()


class CrossGratingEditor(QWidget):
    changed = pyqtSignal()

    def __init__(self, *, title: str = "Cross Grating", allow_overlay: bool = True, parent=None):
        super().__init__(parent)
        self.allow_overlay = allow_overlay
        self._value_store = {}
        self._current_model = ""
        self.param_inputs = {}
        self.child_editors = {}

        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)

        self.group = QGroupBox(title)
        self.group_layout = QFormLayout()
        self.group.setLayout(self.group_layout)
        outer.addWidget(self.group)

        self.model_combo = QComboBox()
        for name in list_cross_grating_models(allow_overlay=allow_overlay):
            self.model_combo.addItem(name)
        self.group_layout.addRow("Model:", self.model_combo)

        self.fields_container = QWidget()
        self.fields_layout = QFormLayout()
        self.fields_container.setLayout(self.fields_layout)
        self.group_layout.addRow(self.fields_container)

        self.model_combo.currentTextChanged.connect(self._on_model_changed)
        initial_model = self.model_combo.currentText() or next(iter(CROSS_GRATING_SPECS))
        self._rebuild_fields(initial_model)
        self.custom_profile_binding = CustomProfileBinding(self, 2, self.to_node, self.from_node,
            lambda: grating_editor_context(self), self._lock_custom_profile_fields, self.changed.emit,
            lambda: grating_editor_media(self))
        self.custom_profile_btn = self.custom_profile_binding.button
        self.group_layout.addRow(self.custom_profile_binding.controls)

    def _lock_custom_profile_fields(self, active):
        for widget in self.param_inputs.values():
            widget.setReadOnly(active)
        for button in self.fields_container.findChildren(QPushButton):
            if button.text() == "Browse...":
                button.setEnabled(not active)

    def closeEvent(self, event):
        self.custom_profile_binding.close()
        super().closeEvent(event)

    def _clear_fields(self):
        while self.fields_layout.count():
            item = self.fields_layout.takeAt(0)
            widget = item.widget()
            if widget is not None:
                widget.deleteLater()
        self.param_inputs = {}
        self.child_editors = {}

    def _add_param_row(self, field, value: str):
        name = field["name"]
        line_edit = QLineEdit()
        line_edit.setText(value)
        line_edit.textChanged.connect(lambda _value: self.changed.emit())
        self.param_inputs[name] = line_edit
        if name == "filename":
            row_widget = QWidget()
            row_layout = QHBoxLayout(row_widget)
            row_layout.setContentsMargins(0, 0, 0, 0)
            row_layout.addWidget(line_edit)
            browse_btn = QPushButton("Browse...")
            browse_btn.clicked.connect(lambda: self._browse_for_file(name))
            row_layout.addWidget(browse_btn)
            self.fields_layout.addRow(f"{field['label']}:", row_widget)
        else:
            self.fields_layout.addRow(f"{field['label']}:", line_edit)

    def _browse_for_file(self, field_name: str):
        fname, _ = QFileDialog.getOpenFileName(self, "Select grating file", "", "All Files (*)")
        if fname:
            widget = self.param_inputs.get(field_name)
            if widget is not None:
                widget.setText(fname)

    def _store_current_state(self):
        if self._current_model:
            self._value_store[self._current_model] = self.to_node()

    def _rebuild_fields(self, model_name: str, node=None):
        choices = list_cross_grating_models(allow_overlay=self.allow_overlay)
        if model_name not in choices:
            model_name = choices[0]
        self._clear_fields()
        state = node if node is not None else self._value_store.get(model_name)
        coerced = coerce_cross_grating(state, allow_overlay=self.allow_overlay, default_model=model_name)
        self._current_model = coerced["model"]
        self._value_store[self._current_model] = coerced
        for field in CROSS_GRATING_SPECS[self._current_model]["fields"]:
            if field["kind"] == "param":
                self._add_param_row(field, coerced["params"][field["name"]])
                continue
            child_name = field["name"]
            child_value = coerced["children"][child_name]
            if field["child_kind"] == "one_d":
                editor = OneDGratingEditor(title=field["label"])
                editor.from_node(child_value)
            else:
                editor = CrossGratingEditor(title=field["label"], allow_overlay=False)
                editor.from_node(child_value)
            self.child_editors[child_name] = editor
            editor.changed.connect(self.changed.emit)
            self.fields_layout.addRow(editor)

    def _on_model_changed(self, model_name: str):
        self.custom_profile_binding.detach()
        self._store_current_state()
        self._rebuild_fields(model_name)
        self.changed.emit()

    def to_node(self, ensure=True):
        if ensure and hasattr(self, "custom_profile_binding"):
            self.custom_profile_binding.ensure()
        choices = list_cross_grating_models(allow_overlay=self.allow_overlay)
        fallback = choices[0] if choices else "OneD_CrossGrating"
        model = self._current_model or self.model_combo.currentText() or fallback
        params = {}
        children = {}
        previous_children = {}
        for field in CROSS_GRATING_SPECS[model]["fields"]:
            if field["kind"] == "param":
                widget = self.param_inputs.get(field["name"])
                value = widget.text().strip() if widget is not None else field["default"]
                params[field["name"]] = value
                continue
            editor = self.child_editors.get(field["name"])
            if editor is not None:
                previous_children[field["name"]] = editor.to_node(ensure=False)
                children[field["name"]] = editor.to_node(ensure=ensure)
        node = {
            "kind": "cross",
            "model": model,
            "params": params,
            "children": previous_children,
        }
        if ensure and children:
            from custom_profile import remap_boundary_media

            node = remap_boundary_media(node, [(previous_children[name], child) for name, child in children.items()])
            for name in ("medium_i", "medium_t"):
                widget = self.param_inputs.get(name)
                if widget is not None and widget.text() != node["params"].get(name, widget.text()):
                    blocked = widget.blockSignals(True)
                    widget.setText(node["params"][name])
                    widget.blockSignals(blocked)
        node["children"] = children
        if ensure and hasattr(self, "custom_profile_binding") and has_editor_profiles(self):
            return normalize_profile_boundaries(node)
        return node

    def from_node(self, node):
        coerced = coerce_cross_grating(node, allow_overlay=self.allow_overlay)
        self._value_store[coerced["model"]] = coerced
        idx = self.model_combo.findText(coerced["model"])
        if idx >= 0:
            previous = self.model_combo.blockSignals(True)
            self.model_combo.setCurrentIndex(idx)
            self.model_combo.blockSignals(previous)
        self._rebuild_fields(coerced["model"], node=coerced)
        self.changed.emit()


class BRDFForm(QWidget):

    DIRECTION_CODES = {
        "Forward Reflection": "0",
        "Forward Transmission": "1",
        "Backward Reflection": "2",
        "Backward Transmission": "3",
    }

    def __init__(self):
        super().__init__()
        self._geometry_viewer = None
        self._cross_profile_state = {}
        self._last_parameter_model = None
        self.custom_profile_binding = CustomProfileBinding(self, 1, self._custom_profile_node,
            self._set_custom_profile_node, self._custom_profile_context, self._lock_custom_profile_fields,
            self._geometry_source_changed, lambda: {"medium_i": "(1,0)", "medium_t": self.substrate.text()})
        self.custom_profile_btn = self.custom_profile_binding.button
        self.custom_profile_binding.controls.hide()

        self.main_layout = QHBoxLayout(self)
        self.setLayout(self.main_layout)


        self.form_layout = QVBoxLayout()
        form_widget = QWidget()
        form_widget.setLayout(self.form_layout)


        self.figure = Figure(figsize=(6, 5), tight_layout=True)
        self.canvas = FigureCanvas(self.figure)
        self.canvas.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)

        plot_layout = QVBoxLayout()
        quantity_row = QHBoxLayout()
        self.plot_quantity = QComboBox()
        try:
            plot_module = importlib.import_module("brdfplot")
        except ImportError:
            plot_module = None
        for key, label in getattr(plot_module, "PLOT_QUANTITIES", (("auto", "Automatic"),)):
            self.plot_quantity.addItem(label, key)
        self.plot_quantity.setToolTip("Choose a quantity from the saved angular result without running the solver again.")
        quantity_row.addWidget(QLabel("Plot:"))
        quantity_row.addWidget(self.plot_quantity, 1)
        plot_layout.addLayout(quantity_row)
        axes_row = QHBoxLayout()
        self.plot_x_axis = QComboBox()
        for label, column in (("Auto", None), ("Scattering angle θs", 0), ("Azimuth φs", 1)):
            self.plot_x_axis.addItem(label, column)
        self.plot_x_axis.setToolTip("These results sweep scattering and azimuth angles at one wavelength.")
        self.plot_y_scale = QComboBox()
        for label, scale in (("Auto", None), ("Linear", False), ("Log", True)):
            self.plot_y_scale.addItem(label, scale)
        self.plot_y_scale.setToolTip("Polarization quantities always use a linear scale so signs and zero remain visible.")
        axes_row.addWidget(QLabel("X angle:"))
        axes_row.addWidget(self.plot_x_axis, 1)
        axes_row.addWidget(QLabel("Y scale:"))
        axes_row.addWidget(self.plot_y_scale)
        plot_layout.addLayout(axes_row)
        self.plot_description = QLabel()
        self.plot_description.setWordWrap(True)
        plot_layout.addWidget(self.plot_description)
        plot_layout.addWidget(self.canvas, 1)
        plot_widget = QWidget()
        plot_widget.setLayout(plot_layout)

        self.main_layout.addWidget(form_widget, 1)
        self.main_layout.addWidget(plot_widget, 1)


        angles_group = QGroupBox("Angular Parameters")
        angle_layout = QFormLayout()
        self.incident_angle = QLineEdit("45")
        self.scatter_start = QLineEdit("45")
        self.scatter_end = QLineEdit("45")
        self.scatter_step = QLineEdit("1")
        self.azimuth_start = QLineEdit("0")
        self.azimuth_end = QLineEdit("180")
        self.azimuth_step = QLineEdit("2")
        angle_layout.addRow("Incident Angle (°):", self.incident_angle)
        angle_layout.addRow("Scattering Start (°):", self.scatter_start)
        angle_layout.addRow("Scattering End (°):", self.scatter_end)
        angle_layout.addRow("Scattering Step (°):", self.scatter_step)
        angle_layout.addRow("Azimuth Start (°):", self.azimuth_start)
        angle_layout.addRow("Azimuth End (°):", self.azimuth_end)
        angle_layout.addRow("Azimuth Step (°):", self.azimuth_step)
        angles_group.setLayout(angle_layout)
        self.form_layout.addWidget(angles_group)


        model_group = QGroupBox("Model Selection")
        model_layout = QFormLayout()
        self.family_selector = QComboBox()
        self.family_selector.addItems([
            "Roughness_BRDF_Model", "Facet_BRDF_Model", "Lambertian_BRDF_Model",
            "Local_BRDF_Model", "Instrument_BRDF_Model", "First_Diffuse_BRDF_Model",
            "Two_Source_BRDF_Model", "Three_Source_BRDF_Model", "Four_Source_BRDF_Model",
            "Transmit_BRDF_Model", "RCW_BRDF_Model", "CrossRCW_BRDF_Model",
            "ZernikeExpansion_BRDF_Model", "Polydisperse_Sphere_BRDF_Model",
        ])
        self.subclass_selector = QComboBox()

        self.subclass_selector.currentTextChanged.connect(self._on_model_changed)
        self.family_selector.currentTextChanged.connect(self.update_subclasses)
        model_layout.addRow("Model Family:", self.family_selector)
        model_layout.addRow("Model Type:", self.subclass_selector)
        model_group.setLayout(model_layout)
        self.form_layout.addWidget(model_group)


        param_group = QGroupBox("Model Parameters")
        param_layout = QFormLayout()
        self.wavelength = QLineEdit("0.532")
        self.substrate = QLineEdit("(4.05,0.05)")
        self.direction = QComboBox()
        self.direction.addItems(list(self.DIRECTION_CODES.keys()))

        param_layout.addRow("Wavelength (µm):", self.wavelength)
        param_layout.addRow("Substrate (n or file):", self.substrate)
        param_layout.addRow("Direction:", self.direction)
        param_group.setLayout(param_layout)
        self.form_layout.addWidget(param_group)


        self.model_params_group = QGroupBox("Selected Model Parameters")
        self.model_params_layout = QFormLayout()
        self.model_params_group.setLayout(self.model_params_layout)
        self.form_layout.addWidget(self.model_params_group)

        self.update_subclasses()
        self.subclass_selector.currentTextChanged.connect(self.populate_model_params)
        self.direction.currentTextChanged.connect(self._on_direction_changed)


        self.psd_group = QGroupBox("PSD Function")
        psd_layout = QVBoxLayout()

        self.psd_function = QComboBox()
        self.psd_function.addItems([
            "Unit_PSD_Function", "ABC_PSD_Function", "Fractal_PSD_Function",
            "Gaussian_PSD_Function", "Elliptical_Mesa_PSD_Function",
            "Rectangular_Mesa_PSD_Function", "Triangular_Mesa_PSD_Function",
            "Rectangular_Pyramid_PSD_Function", "Triangular_Pyramid_PSD_Function",
            "Parabolic_Dimple_PSD_Function",
        ])
        self.psd_function.currentTextChanged.connect(self.update_psd_parameters)
        psd_layout.addWidget(self.psd_function)


        self.psd_param_unit = QWidget()
        self.psd_param_unit.setLayout(QFormLayout())


        self.psd_param_abc = QWidget()
        abc_layout = QFormLayout()
        self.psd_A = QLineEdit("0.01")
        self.psd_B = QLineEdit("362")
        self.psd_C = QLineEdit("2.5")
        abc_layout.addRow("A [µm⁴]:", self.psd_A)
        abc_layout.addRow("B [µm]:", self.psd_B)
        abc_layout.addRow("C [-]:", self.psd_C)
        self.psd_param_abc.setLayout(abc_layout)


        self.psd_param_fractal = QWidget()
        fractal_layout = QFormLayout()
        self.psd_fractalAmp = QLineEdit("0.01")
        self.psd_fractalExp = QLineEdit("2.5")
        fractal_layout.addRow("Amplitude A [µm⁴]:", self.psd_fractalAmp)
        fractal_layout.addRow("Exponent γ [-]:", self.psd_fractalExp)
        self.psd_param_fractal.setLayout(fractal_layout)


        self.psd_param_gaussian = QWidget()
        gaussian_layout = QFormLayout()
        self.psd_sigma = QLineEdit("0.05")
        self.psd_lc = QLineEdit("5")
        gaussian_layout.addRow("Std. Dev. σ [µm]:", self.psd_sigma)
        gaussian_layout.addRow("Correlation Lc [µm]:", self.psd_lc)
        self.psd_param_gaussian.setLayout(gaussian_layout)


        self.psd_param_elliptical_mesa = QWidget()
        el_layout = QFormLayout()
        self.psd_ellipticalX = QLineEdit("1.5")
        self.psd_ellipticalY = QLineEdit("1.5")
        self.psd_mesaHeight = QLineEdit("0.01")
        self.psd_mesaDensity = QLineEdit("0.01")
        el_layout.addRow("X-axis [µm]:", self.psd_ellipticalX)
        el_layout.addRow("Y-axis [µm]:", self.psd_ellipticalY)
        el_layout.addRow("Height [µm]:", self.psd_mesaHeight)
        el_layout.addRow("Density [µm⁻²]:", self.psd_mesaDensity)
        self.psd_param_elliptical_mesa.setLayout(el_layout)


        self.psd_param_rectangular_mesa = QWidget()
        rm_layout = QFormLayout()
        self.psd_rectLenX = QLineEdit("1.5")
        self.psd_rectLenY = QLineEdit("1.5")
        self.psd_rectHeight = QLineEdit("0.01")
        self.psd_rectDensity = QLineEdit("0.01")
        rm_layout.addRow("Length X [µm]:", self.psd_rectLenX)
        rm_layout.addRow("Length Y [µm]:", self.psd_rectLenY)
        rm_layout.addRow("Height [µm]:", self.psd_rectHeight)
        rm_layout.addRow("Density [µm⁻²]:", self.psd_rectDensity)
        self.psd_param_rectangular_mesa.setLayout(rm_layout)


        self.psd_param_triangular_mesa = QWidget()
        tm_layout = QFormLayout()
        self.psd_triSide = QLineEdit("1.5")
        self.psd_triHeight = QLineEdit("0.01")
        self.psd_triDensity = QLineEdit("0.01")
        tm_layout.addRow("Side [µm]:", self.psd_triSide)
        tm_layout.addRow("Height [µm]:", self.psd_triHeight)
        tm_layout.addRow("Density [µm⁻²]:", self.psd_triDensity)
        self.psd_param_triangular_mesa.setLayout(tm_layout)


        self.psd_param_rectangular_pyramid = QWidget()
        rp_layout = QFormLayout()
        self.psd_pyrLenX = QLineEdit("1.5")
        self.psd_pyrLenY = QLineEdit("1.5")
        self.psd_pyrHeight = QLineEdit("0.01")
        self.psd_pyrDensity = QLineEdit("0.01")
        rp_layout.addRow("Base Length X [µm]:", self.psd_pyrLenX)
        rp_layout.addRow("Base Length Y [µm]:", self.psd_pyrLenY)
        rp_layout.addRow("Height [µm]:", self.psd_pyrHeight)
        rp_layout.addRow("Density [µm⁻²]:", self.psd_pyrDensity)
        self.psd_param_rectangular_pyramid.setLayout(rp_layout)


        self.psd_param_triangular_pyramid = QWidget()
        tp_layout = QFormLayout()
        self.psd_tpSide = QLineEdit("1.5")
        self.psd_tpHeight = QLineEdit("0.01")
        self.psd_tpDensity = QLineEdit("0.01")
        tp_layout.addRow("Base Side [µm]:", self.psd_tpSide)
        tp_layout.addRow("Height [µm]:", self.psd_tpHeight)
        tp_layout.addRow("Density [µm⁻²]:", self.psd_tpDensity)
        self.psd_param_triangular_pyramid.setLayout(tp_layout)


        self.psd_param_parabolic_dimple = QWidget()
        pd_layout = QFormLayout()
        self.psd_pdAxisX = QLineEdit("1.5")
        self.psd_pdAxisY = QLineEdit("1.5")
        self.psd_pdHeight = QLineEdit("0.01")
        self.psd_pdDensity = QLineEdit("0.01")
        pd_layout.addRow("Axis X [µm]:", self.psd_pdAxisX)
        pd_layout.addRow("Axis Y [µm]:", self.psd_pdAxisY)
        pd_layout.addRow("Height [µm]:", self.psd_pdHeight)
        pd_layout.addRow("Density [µm⁻²]:", self.psd_pdDensity)
        self.psd_param_parabolic_dimple.setLayout(pd_layout)


        for w in [
            self.psd_param_unit, self.psd_param_abc, self.psd_param_fractal,
            self.psd_param_gaussian, self.psd_param_elliptical_mesa,
            self.psd_param_rectangular_mesa, self.psd_param_triangular_mesa,
            self.psd_param_rectangular_pyramid, self.psd_param_triangular_pyramid,
            self.psd_param_parabolic_dimple,
        ]:
            w.hide()
            psd_layout.addWidget(w)

        self.psd_group.setLayout(psd_layout)
        self.form_layout.addWidget(self.psd_group)


        ctrl_row = QHBoxLayout()
        self.run_btn = QPushButton("Run BRDFProg")
        self.clear_btn = QPushButton("Clear Plot")
        ctrl_row.addWidget(self.run_btn)
        ctrl_row.addWidget(self.clear_btn)
        self.open_output_btn = QPushButton("Open Last Output")
        self.open_input_btn = QPushButton("Open Last Input")
        ctrl_row.addWidget(self.open_output_btn)
        ctrl_row.addWidget(self.open_input_btn)
        self.form_layout.addLayout(ctrl_row)
        self.geometry_btn = QPushButton("View 3D Geometry")
        self.geometry_btn.clicked.connect(self.show_geometry_viewer)
        self.form_layout.addWidget(self.geometry_btn)
        for field in (self.wavelength, self.substrate, self.incident_angle):
            field.textChanged.connect(self._geometry_source_changed)
        self.direction.currentIndexChanged.connect(self._geometry_source_changed)
        self._geometry_source_changed()


        self.output_box = QTextEdit()
        self.output_box.setReadOnly(True)
        self.form_layout.addWidget(QLabel("Log:"))
        self.form_layout.addWidget(self.output_box)


        self.run_btn.clicked.connect(self.run_brdfprog)
        self.clear_btn.clicked.connect(self.clear_plot)
        self.open_output_btn.clicked.connect(self.open_last_output)
        self.open_input_btn.clicked.connect(self.open_last_input)


        self.wavelength.textChanged.connect(
            lambda txt: self._sync_general_to_model("lambda", txt)
        )
        self.substrate.textChanged.connect(
            lambda txt: self._sync_general_to_model("substrate", txt)
        )


        self.update_subclasses()
        self.update_psd_parameters(self.psd_function.currentText())


        self.last_stdout_path = None
        self.last_input_path = None
        self.last_csv_path = None
        self.last_output_meta = None
        self.last_plot_success = False
        self.cross_grating_editor = None
        self._cross_grating_state = build_default_cross_grating()
        for control in (self.plot_quantity, self.plot_x_axis, self.plot_y_scale):
            control.currentIndexChanged.connect(self._replot_last_output)
        self._update_plot_controls()


    def clear_plot(self):
        self.last_plot_success = False
        self.figure.clear()
        self.canvas.draw()

    def _update_plot_controls(self, mod=None, meta=None):
        if mod is None:
            try:
                mod = importlib.import_module("brdfplot")
            except ImportError:
                return
        quantity = self.plot_quantity.currentData() or "auto"
        if quantity == "auto":
            recommend = getattr(mod, "recommended_quantity", None)
            quantity = recommend(meta) if callable(recommend) else "intensity"
        self.plot_y_scale.setEnabled(quantity == "intensity")
        describe = getattr(mod, "quantity_description", None)
        if callable(describe):
            self.plot_description.setText(describe(quantity, meta))
        else:
            self.plot_description.setText("Angular results at the saved wavelength.")

    def _replot_last_output(self, *_):
        self._update_plot_controls(meta=getattr(self, "last_output_meta", None))
        path = getattr(self, "last_csv_path", None)
        if path and Path(path).is_file():
            self.render_with_external(path)


    MODEL_FAMILIES = {

        "Roughness_BRDF_Model": [
            "Microroughness_BRDF_Model",
            "Correlated_Roughness_BRDF_Model",
            "Two_Face_BRDF_Model",
            "Roughness_Stack_BRDF_Model",
            "Correlated_Roughness_Stack_BRDF_Model",
            "Uncorrelated_Roughness_Stack_BRDF_Model",
            "Growth_Roughness_Stack_BRDF_Model",
        ],

        "Facet_BRDF_Model": [
            "Shadowed_Facet_BRDF_Model",
            "Subsurface_Facet_BRDF_Model",
            "Facet_BRDF_Model",
        ],

        "Lambertian_BRDF_Model": [
            "Diffuse_Subsurface_BRDF_Model",
            "Lambertian_BRDF_Model",
        ],

        "Local_BRDF_Model": [
            "Rayleigh_Defect_BRDF_Model",
            "OneLayer_BRDF_Model",
            "Rayleigh_Stack_BRDF_Model",
            "Double_Interaction_BRDF_Model",
            "Subsurface_Particle_BRDF_Model",
            "Bobbert_Vlieger_BRDF_Model",
            "Axisymmetric_Particle_BRDF_Model",
            "Subsurface_Bobbert_Vlieger_BRDF_Model",
            "Subsurface_Axisymmetric_Particle_BRDF_Model",
        ],

        "Instrument_BRDF_Model": [
            "Rayleigh_Instrument_BRDF_Model",
            "Finite_Aperture_Instrument_BRDF_Model",
            "Focussed_Beam_Instrument_BRDF_Model",
        ],
    }

    def update_subclasses(self):







        current_family = self.family_selector.currentText()
        items = self.MODEL_FAMILIES.get(current_family, [])
        self.subclass_selector.clear()
        if items:
            self.subclass_selector.addItems(items)
        else:

            self.subclass_selector.addItems([current_family])


        self.populate_model_params()


    MODEL_PARAM_SPECS = {
        'Microroughness_BRDF_Model': [
            ('lambda', 'double', '0.532'),
            ('substrate', 'dielectric_function', '(4.05,0.05)'),
            ('type', 'int', '0'),
            ('psd', 'PSD_Function_Ptr', 'ABC_PSD_Function'),
        ],
        'Correlated_Roughness_BRDF_Model': [
            ('lambda', 'double', '0.532'),
            ('substrate', 'dielectric_function', '(4.05,0.05)'),
            ('type', 'int', '0'),
            ('psd', 'PSD_Function_Ptr', 'ABC_PSD_Function'),
            ('thickness', 'double', '0.05'),
            ('film', 'dielectric_function', '(1.46,0.05)'),
        ],
        'Two_Face_BRDF_Model': [
            ('lambda', 'double', '0.532'),
            ('substrate', 'dielectric_function', '(4.05,0.05)'),
            ('type', 'int', '0'),
            ('psd', 'PSD_Function_Ptr', 'ABC_PSD_Function'),
            ('film', 'dielectric_function', '(1.46,0)'),
            ('thickness', 'double', '0.05'),
            ('face', 'int', '1'),
        ],
        'Roughness_Stack_BRDF_Model': [
            ('lambda', 'double', '0.532'),
            ('substrate', 'dielectric_function', '(4.05,0.05)'),
            ('type', 'int', '0'),
            ('psd', 'PSD_Function_Ptr', 'ABC_PSD_Function'),
            ('stack', 'StackModel_Ptr', 'No_StackModel'),
            ('this_layer', 'int', '0'),
        ],
        'Correlated_Roughness_Stack_BRDF_Model': [
            ('lambda', 'double', '0.532'),
            ('substrate', 'dielectric_function', '(4.05,0.05)'),
            ('type', 'int', '0'),
            ('psd', 'PSD_Function_Ptr', 'ABC_PSD_Function'),
            ('stack', 'StackModel_Ptr', 'No_StackModel'),
        ],
        'Uncorrelated_Roughness_Stack_BRDF_Model': [
            ('lambda', 'double', '0.532'),
            ('substrate', 'dielectric_function', '(4.05,0.05)'),
            ('type', 'int', '0'),
            ('psd', 'PSD_Function_Ptr', 'ABC_PSD_Function'),
            ('stack', 'StackModel_Ptr', 'No_StackModel'),
        ],
        'Growth_Roughness_Stack_BRDF_Model': [
            ('lambda', 'double', '0.532'),
            ('substrate', 'dielectric_function', '(4.05,0.05)'),
            ('type', 'int', '0'),
            ('psd', 'PSD_Function_Ptr', 'ABC_PSD_Function'),
            ('stack', 'StackModel_Ptr', 'No_StackModel'),
            ('relaxation', 'double', '0.05'),
            ('exponent', 'double', '2'),
            ('intrinsic', 'PSD_Function_Ptr', 'ABC_PSD_Function'),
        ],
        'Shadowed_Facet_BRDF_Model': [
            ('lambda', 'double', '0.532'),
            ('substrate', 'dielectric_function', '(4.05,0.05)'),
            ('type', 'int', '0'),
            ('sdf', 'Slope_Distribution_Function_Ptr', 'Exponential_Slope_Distribution_Function'),
            ('stack', 'StackModel_Ptr', 'No_StackModel'),
            ('shadow', 'Shadow_Function_Ptr', 'Torrance_Sparrow_Shadow_Function'),
        ],
        'Subsurface_Facet_BRDF_Model': [
            ('lambda', 'double', '0.532'),
            ('substrate', 'dielectric_function', '(4.05,0.05)'),
            ('type', 'int', '0'),
            ('sdf', 'Slope_Distribution_Function_Ptr', 'Exponential_Slope_Distribution_Function'),
            ('stack', 'StackModel_Ptr', 'No_StackModel'),
            ('overcoat', 'dielectric_function', '(1.59,0)'),
            ('overcoat_films', 'StackModel_Ptr', 'No_StackModel'),
        ],
        'Lambertian_BRDF_Model': [
            ('lambda', 'double', '0.532'),
            ('substrate', 'dielectric_function', '(4.05,0.05)'),
            ('type', 'int', '0'),
            ('reflectance', 'Reflectance_Ptr', 'Table_Reflectance'),
        ],
        'Diffuse_Subsurface_BRDF_Model': [
            ('lambda', 'double', '0.532'),
            ('substrate', 'dielectric_function', '(4.05,0.05)'),
            ('type', 'int', '0'),
            ('reflectance', 'Reflectance_Ptr', 'Table_Reflectance'),
            ('stack', 'StackModel_Ptr', 'No_StackModel'),
        ],
        'Rayleigh_Defect_BRDF_Model': [
            ('lambda', 'double', '0.532'),
            ('substrate', 'dielectric_function', '(4.05,0.05)'),
            ('type', 'int', '0'),
            ('density', 'double', '1'),
            ('stack', 'StackModel_Ptr', 'No_StackModel'),
            ('radius', 'double', '0.001'),
            ('distance', 'double', '0'),
            ('defect', 'dielectric_function', '(1.0,0.0)'),
        ],
        'OneLayer_BRDF_Model': [
            ('lambda', 'double', '0.532'),
            ('substrate', 'dielectric_function', '(4.05,0.05)'),
            ('type', 'int', '0'),
            ('density', 'double', '1'),
            ('radius', 'double', '0.01'),
            ('defect', 'dielectric_function', '(1.0,0.0)'),
            ('film', 'dielectric_function', '(1.59,0.0)'),
            ('tau', 'double', '0.05'),
            ('depth', 'double', '0.00'),
        ],
        'Rayleigh_Stack_BRDF_Model': [
            ('lambda', 'double', '0.532'),
            ('substrate', 'dielectric_function', '(4.05,0.05)'),
            ('type', 'int', '0'),
            ('density', 'double', '1'),
            ('stack', 'StackModel_Ptr', 'No_StackModel'),
            ('depth', 'double', '0'),
            ('radius', 'double', '0.01'),
            ('sphere', 'dielectric_function', '(1,0)'),
        ],
        'Double_Interaction_BRDF_Model': [
            ('lambda', 'double', '0.532'),
            ('substrate', 'dielectric_function', '(4.05,0.05)'),
            ('type', 'int', '0'),
            ('density', 'double', '1'),
            ('stack', 'StackModel_Ptr', 'No_StackModel'),
            ('distance', 'double', '0.05'),
            ('scatterer', 'Free_Space_Scatterer_Ptr', 'MieScatterer'),
            ('alpha', 'double', '0'),
            ('beta', 'double', '0'),
        ],
        'Subsurface_Particle_BRDF_Model': [
            ('lambda', 'double', '0.532'),
            ('substrate', 'dielectric_function', '(4.05,0.05)'),
            ('type', 'int', '0'),
            ('density', 'double', '1'),
            ('stack', 'StackModel_Ptr', 'No_StackModel'),
            ('depth', 'double', '0'),
            ('scatterer', 'Free_Space_Scatterer_Ptr', 'MieScatterer'),
            ('alpha', 'double', '0'),
            ('beta', 'double', '0'),
        ],
        'Bobbert_Vlieger_BRDF_Model': [
            ('lambda', 'double', '0.532'),
            ('substrate', 'dielectric_function', '(4.05,0.05)'),
            ('type', 'int', '0'),
            ('density', 'double', '1'),
            ('sphere', 'dielectric_function', '(1.59,0)'),
            ('radius', 'double', '0.05'),
            ('spherecoat', 'StackModel_Ptr', 'No_StackModel'),
            ('stack', 'StackModel_Ptr', 'No_StackModel'),
            ('delta', 'double', '0'),
            ('lmax', 'int', '0'),
            ('order', 'int', '-1'),
            ('Norm_Inc_Approx', 'int', '0'),
            ('improve', 'int', '3'),
        ],
        'Axisymmetric_Particle_BRDF_Model': [
            ('lambda', 'double', '0.532'),
            ('substrate', 'dielectric_function', '(4.05,0.05)'),
            ('type', 'int', '0'),
            ('density', 'double', '1'),
            ('Shape', 'Axisymmetric_Shape_Ptr', 'Ellipsoid_Axisymmetric_Shape'),
            ('particle', 'dielectric_function', '(1.59,0)'),
            ('stack', 'StackModel_Ptr', 'No_StackModel'),
            ('delta', 'double', '0'),
            ('lmax', 'int', '0'),
            ('mmax', 'int', '0'),
            ('order', 'int', '-1'),
            ('Norm_Inc_Approx', 'int', '0'),
            ('improve', 'int', '3'),
        ],
        'Subsurface_Bobbert_Vlieger_BRDF_Model': [
            ('lambda', 'double', '0.532'),
            ('substrate', 'dielectric_function', '(4.05,0.05)'),
            ('type', 'int', '0'),
            ('density', 'double', '1'),
            ('sphere', 'dielectric_function', '(1.59,0)'),
            ('radius', 'double', '0.05'),
            ('spherecoat', 'StackModel_Ptr', 'No_StackModel'),
            ('stack', 'StackModel_Ptr', 'No_StackModel'),
            ('delta', 'double', '0'),
            ('lmax', 'int', '0'),
            ('order', 'int', '-1'),
            ('Norm_Inc_Approx', 'int', '0'),
            ('improve', 'int', '3'),
        ],
        'Subsurface_Axisymmetric_Particle_BRDF_Model': [
            ('lambda', 'double', '0.532'),
            ('substrate', 'dielectric_function', '(4.05,0.05)'),
            ('type', 'int', '0'),
            ('density', 'double', '1'),
            ('Shape', 'Axisymmetric_Shape_Ptr', 'Ellipsoid_Axisymmetric_Shape'),
            ('particle', 'dielectric_function', '(1.59,0)'),
            ('stack', 'StackModel_Ptr', 'No_StackModel'),
            ('delta', 'double', '0'),
            ('lmax', 'int', '0'),
            ('mmax', 'int', '0'),
            ('order', 'int', '-1'),
            ('Norm_Inc_Approx', 'int', '0'),
            ('improve', 'int', '3'),
        ],
        'Rayleigh_Instrument_BRDF_Model': [
            ('lambda', 'double', '0.532'),
            ('substrate', 'dielectric_function', '(4.05,0.05)'),
            ('type', 'int', '0'),
            ('stack', 'StackModel_Ptr', 'No_StackModel'),
            ('field_of_view', 'double', '1000'),
            ('air', 'dielectric_function', '(1.0002784,0)'),
            ('number_density', 'double', '2.51E-7'),
        ],
        'Finite_Aperture_Instrument_BRDF_Model': [
            ('lambda', 'double', '0.532'),
            ('substrate', 'dielectric_function', '(4.05,0.05)'),
            ('type', 'int', '0'),
            ('model', 'BRDF_Model_Ptr', 'Microroughness_BRDF_Model'),
            ('alpha', 'double', '0'),
            ('integralmode', 'int', '3'),
        ],
        'Focussed_Beam_Instrument_BRDF_Model': [
            ('lambda', 'double', '0.532'),
            ('substrate', 'dielectric_function', '(4.05,0.05)'),
            ('type', 'int', '0'),
            ('model', 'BRDF_Model_Ptr', 'Microroughness_BRDF_Model'),
            ('alpha', 'double', '0'),
            ('integralmode', 'int', '3'),
            ('focal_point', 'double', '1'),
        ],
        'First_Diffuse_BRDF_Model': [
            ('lambda', 'double', '0.532'),
            ('substrate', 'dielectric_function', '(4.05,0.05)'),
            ('type', 'int', '0'),
            ('depoll', 'double', '0'),
            ('depolc', 'double', '0'),
            ('phase_function', 'Polarized_Phase_Function_Ptr', 'Unpolarized_Phase_Function'),
            ('stack', 'StackModel_Ptr', 'No_StackModel'),
        ],
        'Two_Source_BRDF_Model': [
            ('lambda', 'double', '0.532'),
            ('substrate', 'dielectric_function', '(4.05,0.05)'),
            ('type', 'int', '0'),
            ('factor1', 'double', '1'),
            ('source1', 'BRDF_Model_Ptr', 'Microroughness_BRDF_Model'),
            ('factor2', 'double', '1'),
            ('source2', 'BRDF_Model_Ptr', 'Microroughness_BRDF_Model'),
            ('correlation', 'double', '0'),
        ],
        'Three_Source_BRDF_Model': [
            ('lambda', 'double', '0.532'),
            ('substrate', 'dielectric_function', '(4.05,0.05)'),
            ('type', 'int', '0'),
            ('factor1', 'double', '1'),
            ('source1', 'BRDF_Model_Ptr', 'Microroughness_BRDF_Model'),
            ('factor2', 'double', '1'),
            ('source2', 'BRDF_Model_Ptr', 'Microroughness_BRDF_Model'),
            ('factor3', 'double', '1'),
            ('source3', 'BRDF_Model_Ptr', 'Microroughness_BRDF_Model'),
        ],
        'Four_Source_BRDF_Model': [
            ('lambda', 'double', '0.532'),
            ('substrate', 'dielectric_function', '(4.05,0.05)'),
            ('type', 'int', '0'),
            ('factor1', 'double', '1'),
            ('source1', 'BRDF_Model_Ptr', 'Microroughness_BRDF_Model'),
            ('factor2', 'double', '1'),
            ('source2', 'BRDF_Model_Ptr', 'Microroughness_BRDF_Model'),
            ('factor3', 'double', '1'),
            ('source3', 'BRDF_Model_Ptr', 'Microroughness_BRDF_Model'),
            ('factor4', 'double', '1'),
            ('source4', 'BRDF_Model_Ptr', 'Microroughness_BRDF_Model'),
        ],
        'Transmit_BRDF_Model': [
            ('lambda', 'double', '0.532'),
            ('substrate', 'dielectric_function', '(4.05,0.05)'),
            ('type', 'int', '0'),
            ('model', 'BRDF_Model_Ptr', 'Microroughness_BRDF_Model'),
            ('films', 'StackModel_Ptr', 'No_StackModel'),
        ],
        'RCW_BRDF_Model': [
            ('lambda', 'double', '0.532'),
            ('substrate', 'dielectric_function', '(4.05,0.05)'),
            ('type', 'int', '0'),
            ('alpha', 'double', '0.0175'),
            ('order', 'int', '25'),
            ('grating', 'Grating_Ptr', 'Single_Line_Grating'),
        ],
        'CrossRCW_BRDF_Model': [
            ('lambda', 'double', '0.532'),
            ('substrate', 'dielectric_function', '(4.05,0.05)'),
            ('type', 'int', '0'),
            ('alpha', 'double', '0.01'),
            ('order1', 'int', '10'),
            ('order2', 'int', '10'),
            ('grating', 'CrossGrating_Ptr', 'OneD_CrossGrating'),
        ],
        'ZernikeExpansion_BRDF_Model': [
            ('lambda', 'double', '0.532'),
            ('substrate', 'dielectric_function', '(4.05,0.05)'),
            ('type', 'int', '0'),
            ('coefficientfile', 'string', ''),
            ('scale', 'double', '1'),
        ],
        'Polydisperse_Sphere_BRDF_Model': [
            ('lambda', 'double', '0.532'),
            ('substrate', 'dielectric_function', '(4.05,0.05)'),
            ('type', 'int', '0'),
            ('distribution', 'SurfaceParticleSizeDistribution_Ptr', 'Regular_SurfaceParticleSizeDistribution'),
            ('stack', 'StackModel_Ptr', 'No_StackModel'),
            ('particle', 'dielectric_function', '(1.5,0.0)'),
            ('Dstart', 'double', '0.1'),
            ('Dend', 'double', '100'),
            ('Dstep', 'double', '0.01'),
            ('fractional_coverage', 'double', '0'),
            ('antirainbow', 'double', '0'),
        ],
        'Facet_BRDF_Model': [
            ('lambda', 'double', '0.532'),
            ('substrate', 'dielectric_function', '(4.05,0.05)'),
            ('type', 'int', '0'),
            ('sdf', 'Slope_Distribution_Function_Ptr', 'Exponential_Slope_Distribution_Function'),
            ('stack', 'StackModel_Ptr', 'No_StackModel'),
        ],
    }

    def populate_model_params(self):
        model = self.subclass_selector.currentText() or self.family_selector.currentText()
        previous_model = self._last_parameter_model
        self._last_parameter_model = model
        preserved_expression = None
        if previous_model == model == "RCW_BRDF_Model" and self.custom_profile_binding.document is not None:
            preserved_expression = self._custom_profile_node()
        elif previous_model != model:
            self.custom_profile_binding.detach()
        self.custom_profile_binding.controls.hide()

        existing_editor = getattr(self, "cross_grating_editor", None)
        if existing_editor is not None:
            try:
                self._cross_grating_state = existing_editor.to_node()
                self._cross_profile_state = editor_profiles(existing_editor) if previous_model == model else {}
            except Exception:
                pass
        self.cross_grating_editor = None


        while self.model_params_layout.count():
            item = self.model_params_layout.takeAt(0)
            w = item.widget()
            if w is not None and w is not self.custom_profile_binding.controls:
                w.deleteLater()

        model = self.subclass_selector.currentText() or self.family_selector.currentText()
        specs = self.MODEL_PARAM_SPECS.get(model, [])
        self.param_widgets = {}
        for name, dtype, default in specs:

            if name.lower() in {"psd", "lambda", "substrate", "type"}:
                continue
            if model == "CrossRCW_BRDF_Model" and dtype == "CrossGrating_Ptr" and name == "grating":
                editor = CrossGratingEditor(title="Cross Grating")
                editor.from_node(getattr(self, "_cross_grating_state", build_default_cross_grating()))
                self.cross_grating_editor = editor
                restore_editor_profiles(editor, self._cross_profile_state)
                editor.changed.connect(self._geometry_source_changed)
                self.model_params_layout.addRow(editor)
                continue
            le = QLineEdit(str(default))
            if model in {"RCW_BRDF_Model", "CrossRCW_BRDF_Model"} and name == "alpha":
                le.setToolTip("Detector cone half-angle in radians. This input is separate from the plotted polarization coefficient α = Q/I.")
            if model == "RCW_BRDF_Model" and name == "grating":
                le.textChanged.connect(self._geometry_source_changed)
            self.param_widgets[name] = le
            label = f"{name} ({dtype}):"
            self.model_params_layout.addRow(label, le)
            if model == "RCW_BRDF_Model" and name == "grating":
                self.model_params_layout.addRow(self.custom_profile_binding.controls)
                self.custom_profile_binding.controls.show()
        if preserved_expression is not None:
            self._set_custom_profile_node(preserved_expression)
        self.custom_profile_binding.update_controls()
        self._geometry_source_changed()

    def _custom_profile_context(self):
        return {"wavelength_um": finite_float(self.wavelength.text(), "Wavelength", minimum=0, strict_minimum=True),
                "incident_deg": finite_float(self.incident_angle.text(), "Incident angle", minimum=-90, maximum=90),
                "rotation_deg": 0, "configuration": int(self.DIRECTION_CODES[self.direction.currentText()])}

    def _custom_profile_node(self):
        name, params = self._parse_model_expression(self.param_widgets["grating"].text())
        params.setdefault("medium_t", self.substrate.text())
        return self._grating_expression_node(name, params)

    def _set_custom_profile_node(self, node):
        self.param_widgets["grating"].setText(node_expression(node))

    def _lock_custom_profile_fields(self, active):
        widget = getattr(self, "param_widgets", {}).get("grating")
        if widget is not None:
            widget.setReadOnly(active)

    def _profile_substrate(self):
        active = self._current_model_name() == "RCW_BRDF_Model" and self.custom_profile_binding.document is not None
        active = active or (self._current_model_name() == "CrossRCW_BRDF_Model"
            and self.cross_grating_editor is not None and has_editor_profiles(self.cross_grating_editor))
        if active:
            from custom_profile import native_material_value

            return native_material_value(self.substrate.text())
        return self.substrate.text()

    def _ensure_custom_profiles(self):
        if self._current_model_name() == "RCW_BRDF_Model":
            self.custom_profile_binding.ensure()
        elif self._current_model_name() == "CrossRCW_BRDF_Model" and self.cross_grating_editor is not None:
            self.cross_grating_editor.to_node()

    def geometry_request(self):
        self._ensure_custom_profiles()
        from geometry_types import GeometryRequest

        model = self._current_model_name()
        if model not in {"RCW_BRDF_Model", "CrossRCW_BRDF_Model"}:
            raise ValueError("Select an RCW or CrossRCW BRDF model to preview its geometry.")
        substrate = optical_constant(self._profile_substrate(), "Substrate")
        if model == "CrossRCW_BRDF_Model":
            if self.cross_grating_editor is None:
                raise ValueError("Cross grating inputs are unavailable.")
            tree = self.cross_grating_editor.to_node()
            errors = validate_cross_grating(tree)
            if optical_constant(tree["params"]["medium_i"], "Grating incident medium") != "(1,0)":
                raise ValueError("CrossRCW BRDF requires grating medium_i to be vacuum (1,0).")
        else:
            name, params = self._parse_model_expression(self.param_widgets["grating"].text())
            params.setdefault("medium_t", substrate)
            tree = self._grating_expression_node(name, params)
            errors = validate_one_d_grating(tree)
        if errors:
            raise ValueError("Grating validation failed: " + "; ".join(errors))
        if optical_constant(tree["params"]["medium_t"], "Grating transmitted medium") != substrate:
            raise ValueError("Grating medium_t must match the BRDF substrate.")
        wavelength = finite_float(self.wavelength.text(), "Wavelength", minimum=0, strict_minimum=True)
        incident = finite_float(self.incident_angle.text(), "Incident angle", minimum=-90, maximum=90)
        if abs(incident) == 90:
            raise ValueError("Incident angle must be strictly between -90 and 90 degrees.")
        return GeometryRequest(
            grating=tree, wavelength_um=wavelength, incident_deg=incident,
            configuration=int(self.DIRECTION_CODES[self.direction.currentText()]),
            source_label="BRDF · " + ("CrossRCW" if model == "CrossRCW_BRDF_Model" else "RCW"),
        )

    def show_geometry_viewer(self):
        if self._geometry_viewer is None:
            from geometry_viewer import GeometryViewer
            self._geometry_viewer = GeometryViewer(parent=self, request_provider=self.geometry_request)
        self._geometry_viewer.show()
        self._geometry_viewer.raise_()
        self._geometry_viewer.activateWindow()

    def _geometry_source_changed(self, *_):
        if hasattr(self, "geometry_btn"):
            self.geometry_btn.setEnabled(self._current_model_name() in {"RCW_BRDF_Model", "CrossRCW_BRDF_Model"})
        if self._geometry_viewer is not None:
            self._geometry_viewer.source_changed()

    def closeEvent(self, event):
        self.custom_profile_binding.close()
        if self.cross_grating_editor is not None:
            for editor in [self.cross_grating_editor] + self.cross_grating_editor.findChildren(QWidget):
                binding = getattr(editor, "custom_profile_binding", None)
                if binding is not None:
                    binding.close()
        if self._geometry_viewer is not None:
            self._geometry_viewer.close()
        super().closeEvent(event)

    def update_psd_parameters(self, current: str):
        mapping = {
            "Unit_PSD_Function": self.psd_param_unit,
            "ABC_PSD_Function": self.psd_param_abc,
            "Fractal_PSD_Function": self.psd_param_fractal,
            "Gaussian_PSD_Function": self.psd_param_gaussian,
            "Elliptical_Mesa_PSD_Function": self.psd_param_elliptical_mesa,
            "Rectangular_Mesa_PSD_Function": self.psd_param_rectangular_mesa,
            "Triangular_Mesa_PSD_Function": self.psd_param_triangular_mesa,
            "Rectangular_Pyramid_PSD_Function": self.psd_param_rectangular_pyramid,
            "Triangular_Pyramid_PSD_Function": self.psd_param_triangular_pyramid,
            "Parabolic_Dimple_PSD_Function": self.psd_param_parabolic_dimple,
        }
        for w in mapping.values():
            w.hide()
        widget = mapping.get(current)
        if widget is not None:
            widget.show()

    def _sync_general_to_model(self, name: str, text: str):
        widget = getattr(self, "param_widgets", {}).get(name)
        if widget is not None and widget.text() != text:
            widget.setText(text)

    def _on_direction_changed(self, selection: str):
        value = self.DIRECTION_CODES.get(selection)
        if value is None:
            return
        widget = getattr(self, "param_widgets", {}).get("type")
        if widget is not None and widget.text() != value:
            widget.setText(value)

    def _current_model_name(self) -> str:
        return self.subclass_selector.currentText() or self.family_selector.currentText()

    def _current_model_specs(self):
        return self.MODEL_PARAM_SPECS.get(self._current_model_name(), [])

    def _collect_model_params(self) -> dict:
        model = self._current_model_name()
        specs = self.MODEL_PARAM_SPECS.get(model, [])
        params = {}
        for name, dtype, default in specs:
            lowered = name.lower()
            if lowered in {"psd", "lambda", "substrate", "type"}:
                continue
            if model == "CrossRCW_BRDF_Model" and dtype == "CrossGrating_Ptr" and name == "grating":
                continue
            widget = self.param_widgets.get(name)
            if widget is None:
                continue
            text = widget.text().strip()
            params[name] = text
        return params

    def _current_cross_grating_tree(self):
        editor = getattr(self, "cross_grating_editor", None)
        if editor is not None:
            tree = editor.to_node()
            self._cross_grating_state = tree
            return tree
        return getattr(self, "_cross_grating_state", build_default_cross_grating())

    def _append_cross_rcw_input_lines(self, input_lines, *, model_params=None, grating_tree=None):
        params = dict(model_params or {})
        alpha = str(params.get("alpha", "0.01") or "0.01")
        order1 = str(params.get("order1", "10") or "10")
        order2 = str(params.get("order2", "10") or "10")
        tree = grating_tree if grating_tree is not None else self._current_cross_grating_tree()
        errors = validate_cross_grating(tree)
        if errors:
            raise ValueError("Cross grating validation failed:\n- " + "\n- ".join(errors))
        input_lines.extend([alpha, order1, order2])
        input_lines.extend(serialize_cross_grating(tree))
        self._cross_grating_state = tree


    def render_with_external(self, csv_path: str, *, x_col: int = None, y_col: int = None, semilogy: bool = None):
        csv_path = str(csv_path)
        self.last_csv_path = csv_path
        self.last_plot_success = False
        self.figure.clear()
        ax = self.figure.add_subplot(111)

        here = Path(__file__).resolve().parent
        cwd = Path(os.getcwd())
        csv_dir = Path(csv_path).resolve().parent if csv_path else None

        meta = getattr(self, "last_output_meta", None)
        if not meta or meta.get("csv_path") != csv_path:
            meta = self._load_output_meta(csv_path)
        self.last_output_meta = meta

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
                    spec = importlib.util.spec_from_file_location("brdfplot", str(path))
                    if spec and spec.loader:
                        mod = importlib.util.module_from_spec(spec)
                        sys.modules["brdfplot"] = mod
                        spec.loader.exec_module(mod)
                        how = f"spec_from_file_location({path})"
                        return True
            except Exception as e:
                tried.append(f"load {path}: {e}")
            return False

        if not _try_import("brdfplot"):
            if not _try_path(here / "brdfplot.py"):
                if not _try_path(cwd / "brdfplot.py") and csv_dir:
                    _try_path(csv_dir / "brdfplot.py")

        if mod is None:
            self.output_box.append("Could not import brdfplot.py: " + " | ".join(tried))
            self.canvas.draw()
            return
        else:
            self.output_box.append(f"brdfplot resolved via: {how}")

        fn = getattr(mod, "plot_csv", None)
        if not callable(fn):
            self.output_box.append("brdfplot.py found, but it must define plot_csv(ax, csv_path, ...).")
            self.canvas.draw()
            return

        try:
            self._update_plot_controls(mod, meta)
            if x_col is None:
                x_col = self.plot_x_axis.currentData()
            if semilogy is None:
                semilogy = self.plot_y_scale.currentData() if self.plot_y_scale.isEnabled() else False
            kwargs = {"x_col": x_col, "y_col": y_col, "semilogy": semilogy,
                      "meta": meta, "quantity": self.plot_quantity.currentData() or "auto"}
            try:
                parameters = inspect.signature(fn).parameters
            except (TypeError, ValueError):
                parameters = None
            if parameters is not None and not any(value.kind == inspect.Parameter.VAR_KEYWORD
                                                  for value in parameters.values()):
                kwargs = {name: value for name, value in kwargs.items() if name in parameters}
            fn(ax, csv_path, **kwargs)
            self.canvas.draw()
            self.last_plot_success = True
            self.output_box.append("Plot updated via brdfplot.plot_csv")
        except Exception as e:
            self.output_box.append(f"brdfplot render error: {e}")
            self.canvas.draw()

    def _split_top_level_args(self, text: str):
        parts = []
        token = []
        depth = 0
        for ch in text:
            if ch == "(":
                depth += 1
                token.append(ch)
            elif ch == ")":
                if depth == 0:
                    raise ValueError("Unbalanced parentheses in model expression.")
                depth -= 1
                token.append(ch)
            elif ch == "," and depth == 0:
                item = "".join(token).strip()
                if item:
                    parts.append(item)
                token = []
            else:
                token.append(ch)
        if depth:
            raise ValueError("Unbalanced parentheses in model expression.")
        item = "".join(token).strip()
        if item:
            parts.append(item)
        return parts

    def _parse_model_expression(self, text: str):
        raw = (text or "").strip()
        if not raw:
            return "", {}

        m = re.match(r"^([A-Za-z_][A-Za-z0-9_]*)\s*(?:\((.*)\))?\s*$", raw)
        if not m:
            raise ValueError(f"Invalid model expression: {raw}")

        model_name = m.group(1).strip()
        inner = m.group(2)
        if inner is None:
            return model_name, {}

        kwargs = {}
        for item in self._split_top_level_args(inner):
            if "=" not in item:
                raise ValueError(f"Expected name=value in model expression: {item}")
            key, value = item.split("=", 1)
            key, value = key.strip().lower(), value.strip()
            if not key or not value or key in kwargs:
                raise ValueError(f"Invalid or duplicate model parameter: {key}")
            kwargs[key] = value
        return model_name, kwargs

    def _append_scatterer_input_lines(self, input_lines, value: str, model_name: str):
        scatterer, kwargs = self._parse_model_expression(value)
        kwargs.setdefault("lambda", self.wavelength.text())
        if model_name == "Subsurface_Particle_BRDF_Model":
            kwargs.setdefault("medium", self.substrate.text())
        input_lines.append(scatterer)
        self._append_model_parameters(input_lines, scatterer, kwargs)

    def _append_psd_input_lines(self, input_lines):
        current = self.psd_function.currentText()
        input_lines.append(current)
        if current == "Unit_PSD_Function":
            return
        if current == "ABC_PSD_Function":
            input_lines += [self.psd_A.text(), self.psd_B.text(), self.psd_C.text()]
            return
        if current == "Fractal_PSD_Function":
            input_lines += [self.psd_fractalAmp.text(), self.psd_fractalExp.text()]
            return
        if current == "Gaussian_PSD_Function":
            input_lines += [self.psd_sigma.text(), self.psd_lc.text()]
            return
        if current == "Elliptical_Mesa_PSD_Function":
            input_lines += [
                self.psd_ellipticalX.text(),
                self.psd_ellipticalY.text(),
                self.psd_mesaHeight.text(),
                self.psd_mesaDensity.text(),
            ]
            return
        if current == "Rectangular_Mesa_PSD_Function":
            input_lines += [
                self.psd_rectLenX.text(),
                self.psd_rectLenY.text(),
                self.psd_rectHeight.text(),
                self.psd_rectDensity.text(),
            ]
            return
        if current == "Triangular_Mesa_PSD_Function":
            input_lines += [self.psd_triSide.text(), self.psd_triHeight.text(), self.psd_triDensity.text()]
            return
        if current == "Rectangular_Pyramid_PSD_Function":
            input_lines += [
                self.psd_pyrLenX.text(),
                self.psd_pyrLenY.text(),
                self.psd_pyrHeight.text(),
                self.psd_pyrDensity.text(),
            ]
            return
        if current == "Triangular_Pyramid_PSD_Function":
            input_lines += [self.psd_tpSide.text(), self.psd_tpHeight.text(), self.psd_tpDensity.text()]
            return
        if current == "Parabolic_Dimple_PSD_Function":
            input_lines += [
                self.psd_pdAxisX.text(),
                self.psd_pdAxisY.text(),
                self.psd_pdHeight.text(),
                self.psd_pdDensity.text(),
            ]
            return
        raise ValueError(f"Unsupported PSD function: {current}")

    OUTPUT_COLUMNS = ["theta_s", "phi_s", "eta", "DOLP", "DOCP", "DOP", "I"]

    @staticmethod
    def _validated_parameter(name, dtype, value):
        label = f"Model parameter {name}"
        if dtype == "double":
            positive = name in {"lambda", "radius", "Dstart", "Dend", "Dstep"}
            nonnegative = name in {"density", "thickness", "depth", "tau", "number_density", "field_of_view"}
            return str(finite_float(value, label, minimum=0 if positive or nonnegative else None,
                                    strict_minimum=positive))
        if dtype == "int":
            bounds = {"type": (0, 3), "face": (0, 1), "order1": (0, None),
                      "order2": (0, None), "integralmode": (1, 7), "Norm_Inc_Approx": (0, 1)}
            lo, hi = bounds.get(name, (None, None))
            return str(integer(value, label, minimum=lo, maximum=hi))
        if dtype == "dielectric_function":
            return optical_constant(value, label)
        if dtype == "dielectric_stack":
            return ":" + solver_file(single_line(value, label).removeprefix(":"), label)
        if dtype == "string" and name == "coefficientfile":
            return solver_file(value, label)
        return single_line(value, label)

    def _build_input_lines(self):
        self._ensure_custom_profiles()
        angle_specs = [("incident_angle", "Incident angle", -90, 90),
                       ("scatter_start", "Scattering start", -90, 90),
                       ("scatter_end", "Scattering end", -90, 90),
                       ("scatter_step", "Scattering step", 0, None),
                       ("azimuth_start", "Azimuth start", None, None),
                       ("azimuth_end", "Azimuth end", None, None),
                       ("azimuth_step", "Azimuth step", 0, None)]
        angles = {}
        for attr, label, lo, hi in angle_specs:
            angles[attr] = finite_float(getattr(self, attr).text(), label, minimum=lo, maximum=hi,
                                        strict_minimum=attr.endswith("step"))
        for prefix in ("scatter", "azimuth"):
            start, end, step = (angles[f"{prefix}_{name}"] for name in ("start", "end", "step"))
            if end < start:
                raise ValueError(f"{prefix.capitalize()} end must be at least its start.")
            for bound in (start, end):
                radians = math.radians(bound)
                if radians + math.radians(step) <= radians:
                    raise ValueError(f"{prefix.capitalize()} step is too small to advance the scan.")
        model = self._current_model_name()
        params = self._collect_model_params()
        params.update({"lambda": self.wavelength.text(), "substrate": self._profile_substrate(),
                       "type": self.DIRECTION_CODES[self.direction.currentText()]})
        lines = [str(angles[attr]) for attr, *_ in angle_specs]
        lines.append(model)
        self._append_model_parameters(lines, model, params, use_selected_psd=True)
        return lines

    NESTED_MODEL_SPECS = {
        'ABC_PSD_Function': [
            ('A', 'double', '0.01'),
            ('B', 'double', '362'),
            ('C', 'double', '2.5'),
        ],
        'Anisotropic_Exponential_Slope_Distribution_Function': [
            ('sx', 'double', '0.1'),
            ('sy', 'double', '0.01'),
        ],
        'Bimodal_Distribution': [
            ('distributionA', 'Distribution_Ptr', 'Log_Normal_Distribution'),
            ('distributionB', 'Distribution_Ptr', 'Log_Normal_Distribution'),
            ('fractionA', 'double', '0.5'),
        ],
        'CC1246E_SurfaceParticleSizeDistribution': [
            ('cleanliness', 'double', '1000'),
            ('slope', 'double', '0.926'),
        ],
        'Chebyshev_Axisymmetric_Shape': [
            ('npoints', 'int', '100'),
            ('vertical', 'double', '0.05'),
            ('horizontal', 'double', '0.05'),
            ('amplitude', 'double', '0.01'),
            ('n', 'int', '8'),
            ('start', 'double', '0'),
            ('end', 'double', '3.14159265'),
            ('renorm', 'int', '0'),
        ],
        'CoatedMieScatterer': [
            ('lambda', 'double', '0.532'),
            ('medium', 'dielectric_function', '(1,0)'),
            ('radius', 'double', '0.05'),
            ('sphere', 'dielectric_function', '(1.59,0)'),
            ('coating', 'dielectric_function', '(1,0)'),
            ('thickness', 'double', '0'),
            ('nmax', 'int', '0'),
        ],
        'CodeWhitney_Phase_Function': [
            ('phase_function', 'Phase_Function_Ptr', 'Henyey_Greenstein_Phase_Function'),
            ('log_absorption', 'double', '-20'),
            ('log_scattering', 'double', '0'),
            ('plmax', 'double', '1'),
            ('pcmax', 'double', '0'),
            ('s', 'double', '0'),
        ],
        'Conical_Axisymmetric_Shape': [
            ('npoints', 'int', '100'),
            ('radius', 'double', '0.05'),
            ('height', 'double', '0.1'),
            ('offset', 'double', '0'),
            ('corner_base', 'double', '0'),
            ('corner_apex', 'double', '0'),
            ('updown', 'int', '0'),
            ('renorm', 'int', '0'),
        ],
        'Constrained_Phase_Function': [
            ('phase_function', 'Phase_Function_Ptr', 'Henyey_Greenstein_Phase_Function'),
            ('log_absorption', 'double', '-20'),
            ('log_scattering', 'double', '0'),
            ('m11_0', 'double', '1'),
            ('m11_60', 'double', '0'),
            ('m11_120', 'double', '0'),
            ('m11_180', 'double', '0'),
            ('m22_0', 'double', '1'),
            ('m22_60', 'double', '0'),
            ('m22_120', 'double', '0'),
            ('m22_180', 'double', '0'),
            ('m33_0', 'double', '1'),
            ('m33_60', 'double', '0'),
            ('m33_120', 'double', '0'),
            ('m33_180', 'double', '0'),
            ('m10_60', 'double', '0'),
            ('m10_120', 'double', '0'),
            ('m23_60', 'double', '0'),
            ('m23_120', 'double', '0'),
        ],
        'Cylinder_Axisymmetric_Shape': [
            ('npoints', 'int', '100'),
            ('radius', 'double', '0.05'),
            ('length', 'double', '0.1'),
            ('corner', 'double', '0'),
            ('renorm', 'int', '0'),
        ],
        'Distributed_Mie_Phase_Function': [
            ('lambda', 'double', '0.532'),
            ('medium', 'dielectric_function', '(1.33,0)'),
            ('sphere', 'dielectric_function', '(1,0)'),
            ('distribution', 'VolumeParticleSizeDistribution_Ptr', 'Regular_VolumeParticleSizeDistribution'),
            ('integralStart', 'double', '0.1'),
            ('integralEnd', 'double', '100'),
            ('integralStep', 'double', '0.01'),
        ],
        'DoubleFilm_StackModel': [
            ('material1', 'dielectric_function', '(1.5,0)'),
            ('thickness1', 'double', '0.1'),
            ('material2', 'dielectric_function', '(1.5,0)'),
            ('thickness2', 'double', '0.1'),
        ],
        'Double_Conical_Axisymmetric_Shape': [
            ('npoints', 'int', '100'),
            ('height', 'double', '0.1'),
            ('radius', 'double', '0.05'),
            ('corner_apex', 'double', '0'),
            ('corner_waist', 'double', '0'),
            ('renorm', 'int', '0'),
        ],
        'Double_Henyey_Greenstein_Phase_Function': [
            ('g1', 'double', '0.01'),
            ('g2', 'double', '0.01'),
            ('c', 'double', '0.1'),
        ],
        'Double_PSD_Function': [
            ('psd1', 'PSD_Function_Ptr', 'ABC_PSD_Function'),
            ('psd2', 'PSD_Function_Ptr', 'ABC_PSD_Function'),
        ],
        'Double_Reynolds_McCormick_Phase_Function': [
            ('g', 'double', '0.01'),
            ('alpha', 'double', '0.5'),
            ('c', 'double', '0.1'),
        ],
        'Ellipsoid_Axisymmetric_Shape': [
            ('npoints', 'int', '100'),
            ('vertical', 'double', '0.05'),
            ('horizontal', 'double', '0.05'),
            ('offset', 'double', '0'),
        ],
        'Ellipsoid_Slope_Distribution_Function': [
            ('rx', 'double', '10'),
            ('ry', 'double', '10'),
            ('rz', 'double', '1'),
            ('density', 'double', '0.0001'),
        ],
        'Elliptical_Mesa_PSD_Function': [
            ('axisx', 'double', '1.5'),
            ('axisy', 'double', '1.5'),
            ('height', 'double', '0.01'),
            ('density', 'double', '0.01'),
        ],
        'Equation_Reflectance': [
            ('expression', 'string', '1'),
        ],
        'Exponential_Angle_Distribution_Function': [
            ('sigma', 'double', '5'),
        ],
        'Exponential_Slope_Distribution_Function': [
            ('s', 'double', '0.1'),
        ],
        'Fractal_PSD_Function': [
            ('A', 'double', '0.01'),
            ('exponent', 'double', '2.5'),
        ],
        'Gaussian_Angle_Distribution_Function': [
            ('sigma', 'double', '5'),
        ],
        'Gaussian_PSD_Function': [
            ('sigma', 'double', '0.05'),
            ('length', 'double', '5'),
        ],
        'Gaussian_Phase_Function': [
            ('width', 'double', '20'),
            ('background', 'double', '0.1'),
        ],
        'Gaussian_Slope_Distribution_Function': [
            ('s', 'double', '0.1'),
        ],
        'GradedFilm_StackModel': [
            ('start', 'dielectric_function', '(1.5,0)'),
            ('end', 'dielectric_function', '(1,0)'),
            ('thickness', 'double', '0.1'),
            ('steps', 'int', '5'),
        ],
        'Henyey_Greenstein_Phase_Function': [
            ('g', 'double', '0.01'),
        ],
        'Indented_Ellipsoid_Axisymmetric_Shape': [
            ('npoints', 'int', '100'),
            ('vertical', 'double', '0.05'),
            ('horizontal', 'double', '0'),
            ('indent', 'double', '0.01'),
            ('updown', 'int', '0'),
            ('renorm', 'int', '0'),
        ],
        'Isotropic_Phase_Function': [
        ],
        'Kubelka_Munk_Reflectance': [
            ('thickness', 'double', '1'),
            ('absorption', 'Table', '1'),
            ('scattering', 'Table', '1'),
            ('substrate', 'Reflectance_Ptr', 'Table_Reflectance'),
        ],
        'Legendre_Phase_Function': [
            ('c0', 'double', '1'),
            ('c1', 'double', '0.65'),
            ('c2', 'double', '0.42'),
            ('c3', 'double', '0'),
            ('c4', 'double', '0'),
            ('c5', 'double', '0'),
        ],
        'Log_Normal_Distribution': [
            ('sigma', 'double', '1'),
            ('median', 'double', '10'),
        ],
        'Maxwell_Beard_Shadow_Function': [
            ('Omega', 'double', '0.7'),
            ('tau', 'double', '0.25'),
        ],
        'MieScatterer': [
            ('lambda', 'double', '0.532'),
            ('medium', 'dielectric_function', '(1,0)'),
            ('radius', 'double', '0.05'),
            ('sphere', 'dielectric_function', '(1.59,0)'),
        ],
        'Modified_Gamma_Distribution': [
            ('mu', 'double', '8'),
            ('Lambda', 'double', '0.0417'),
            ('gamma', 'double', '3'),
        ],
        'MultilayerCoatedMieScatterer': [
            ('lambda', 'double', '0.532'),
            ('medium', 'dielectric_function', '(1,0)'),
            ('radius', 'double', '0.05'),
            ('sphere', 'dielectric_function', '(1.59,0)'),
            ('stack', 'StackModel_Ptr', 'No_StackModel'),
        ],
        'No_StackModel': [
        ],
        'Normal_Distribution': [
            ('mean', 'double', '10'),
            ('sigma', 'double', '5'),
        ],
        'Optically_Active_Sphere_Scatterer': [
            ('lambda', 'double', '0.532'),
            ('medium', 'dielectric_function', '(1,0)'),
            ('materialleft', 'dielectric_function', '(1.59,0)'),
            ('materialright', 'dielectric_function', '(1.59,0)'),
            ('radius', 'double', '1'),
        ],
        'Parabolic_Dimple_PSD_Function': [
            ('axisx', 'double', '1.5'),
            ('axisy', 'double', '1.5'),
            ('height', 'double', '0.01'),
            ('density', 'double', '0.01'),
        ],
        'RayleighGansSphereScatterer': [
            ('lambda', 'double', '0.532'),
            ('medium', 'dielectric_function', '(1,0)'),
            ('radius', 'double', '0.05'),
            ('sphere', 'dielectric_function', '(1.59,0)'),
        ],
        'RayleighScatterer': [
            ('lambda', 'double', '0.532'),
            ('medium', 'dielectric_function', '(1,0)'),
            ('radius', 'double', '0.05'),
            ('sphere', 'dielectric_function', '(1.59,0)'),
        ],
        'Rayleigh_Phase_Function': [
        ],
        'Rectangular_Mesa_PSD_Function': [
            ('lengthx', 'double', '1.5'),
            ('lengthy', 'double', '1.5'),
            ('height', 'double', '0.01'),
            ('density', 'double', '0.01'),
        ],
        'Rectangular_Pyramid_PSD_Function': [
            ('lengthx', 'double', '1.5'),
            ('lengthy', 'double', '1.5'),
            ('height', 'double', '0.01'),
            ('density', 'double', '0.01'),
        ],
        'Regular_SurfaceParticleSizeDistribution': [
            ('distribution', 'Distribution_Ptr', 'Log_Normal_Distribution'),
            ('numberdensity', 'double', '0.001'),
        ],
        'Regular_VolumeParticleSizeDistribution': [
            ('distribution', 'Distribution_Ptr', 'Log_Normal_Distribution'),
            ('numberdensity', 'double', '0.00001'),
        ],
        'Reynolds_McCormick_Phase_Function': [
            ('g', 'double', '0.01'),
            ('alpha', 'double', '0.5'),
        ],
        'SingleFilm_StackModel': [
            ('material', 'dielectric_function', '(1.5,0)'),
            ('thickness', 'double', '0.1'),
        ],
        'Smith_Shadow_Function': [
            ('w', 'double', '0.35'),
        ],
        'SphericalScatterer_Phase_Function': [
            ('scatterer', 'SphericalScatterer_Ptr', 'MieScatterer'),
            ('volume_density', 'double', '1E-4'),
        ],
        'Stack_StackModel': [
            ('stack', 'dielectric_stack', ''),
        ],
        'TMatrix_Axisymmetric_Scatterer': [
            ('lambda', 'double', '0.532'),
            ('medium', 'dielectric_function', '(1,0)'),
            ('Shape', 'Axisymmetric_Shape_Ptr', 'Ellipsoid_Axisymmetric_Shape'),
            ('particle', 'dielectric_function', '(1.59,0)'),
            ('lmax', 'int', '0'),
            ('mmax', 'int', '0'),
        ],
        'Table_Angle_Distribution_Function': [
            ('T', 'Table', '1'),
        ],
        'Table_Axisymmetric_Shape': [
            ('npoints', 'int', '100'),
            ('table', 'Table', '0.05'),
            ('scale', 'double', '1'),
        ],
        'Table_PSD_Function': [
            ('T', 'Table', '1'),
        ],
        'Table_Phase_Function': [
            ('table', 'Table', '0.07957747'),
        ],
        'Table_Reflectance': [
            ('table', 'Table', '1'),
        ],
        'Table_Shadow_Function': [
            ('T', 'Table', '1'),
        ],
        'Table_Slope_Distribution_Function': [
            ('T', 'Table', '1'),
        ],
        'Torrance_Sparrow_Shadow_Function': [
        ],
        'Triangular_Mesa_PSD_Function': [
            ('side', 'double', '1.5'),
            ('height', 'double', '0.01'),
            ('density', 'double', '0.01'),
        ],
        'Triangular_Pyramid_PSD_Function': [
            ('length', 'double', '1.5'),
            ('height', 'double', '0.01'),
            ('density', 'double', '0.01'),
        ],
        'Two_Slope_Distribution_Function': [
            ('s1', 'Slope_Distribution_Function_Ptr', 'Exponential_Slope_Distribution_Function'),
            ('s2', 'Slope_Distribution_Function_Ptr', 'Exponential_Slope_Distribution_Function'),
            ('fract', 'double', '0.5'),
        ],
        'Unit_PSD_Function': [
        ],
        'Unit_Shadow_Function': [
        ],
        'Unit_Slope_Distribution_Function': [
        ],
        'Unpolarized_Phase_Function': [
            ('phase_function', 'Phase_Function_Ptr', 'Henyey_Greenstein_Phase_Function'),
            ('log_absorption', 'double', '-20'),
            ('log_scattering', 'double', '0'),
        ],
        'Weibull_Distribution': [
            ('scale', 'double', '10'),
            ('exponent', 'double', '4'),
        ],
    }

    NESTED_MODEL_TYPES = {
        'StackModel_Ptr': {
            'DoubleFilm_StackModel',
            'GradedFilm_StackModel',
            'No_StackModel',
            'SingleFilm_StackModel',
            'Stack_StackModel',
        },
        'PSD_Function_Ptr': {
            'ABC_PSD_Function',
            'Double_PSD_Function',
            'Elliptical_Mesa_PSD_Function',
            'Fractal_PSD_Function',
            'Gaussian_PSD_Function',
            'Parabolic_Dimple_PSD_Function',
            'Rectangular_Mesa_PSD_Function',
            'Rectangular_Pyramid_PSD_Function',
            'Table_PSD_Function',
            'Triangular_Mesa_PSD_Function',
            'Triangular_Pyramid_PSD_Function',
            'Unit_PSD_Function',
        },
        'Slope_Distribution_Function_Ptr': {
            'Anisotropic_Exponential_Slope_Distribution_Function',
            'Ellipsoid_Slope_Distribution_Function',
            'Exponential_Angle_Distribution_Function',
            'Exponential_Slope_Distribution_Function',
            'Gaussian_Angle_Distribution_Function',
            'Gaussian_Slope_Distribution_Function',
            'Table_Angle_Distribution_Function',
            'Table_Slope_Distribution_Function',
            'Two_Slope_Distribution_Function',
            'Unit_Slope_Distribution_Function',
        },
        'Shadow_Function_Ptr': {
            'Maxwell_Beard_Shadow_Function',
            'Smith_Shadow_Function',
            'Table_Shadow_Function',
            'Torrance_Sparrow_Shadow_Function',
            'Unit_Shadow_Function',
        },
        'Reflectance_Ptr': {
            'Equation_Reflectance',
            'Kubelka_Munk_Reflectance',
            'Table_Reflectance',
        },
        'Axisymmetric_Shape_Ptr': {
            'Chebyshev_Axisymmetric_Shape',
            'Conical_Axisymmetric_Shape',
            'Cylinder_Axisymmetric_Shape',
            'Double_Conical_Axisymmetric_Shape',
            'Ellipsoid_Axisymmetric_Shape',
            'Indented_Ellipsoid_Axisymmetric_Shape',
            'Table_Axisymmetric_Shape',
        },
        'Free_Space_Scatterer_Ptr': {
            'CoatedMieScatterer',
            'MieScatterer',
            'MultilayerCoatedMieScatterer',
            'Optically_Active_Sphere_Scatterer',
            'RayleighGansSphereScatterer',
            'RayleighScatterer',
            'TMatrix_Axisymmetric_Scatterer',
        },
        'SphericalScatterer_Ptr': {
            'CoatedMieScatterer',
            'MieScatterer',
            'MultilayerCoatedMieScatterer',
            'RayleighGansSphereScatterer',
            'RayleighScatterer',
        },
        'Polarized_Phase_Function_Ptr': {
            'CodeWhitney_Phase_Function',
            'Constrained_Phase_Function',
            'Distributed_Mie_Phase_Function',
            'SphericalScatterer_Phase_Function',
            'Unpolarized_Phase_Function',
        },
        'Phase_Function_Ptr': {
            'Double_Henyey_Greenstein_Phase_Function',
            'Double_Reynolds_McCormick_Phase_Function',
            'Gaussian_Phase_Function',
            'Henyey_Greenstein_Phase_Function',
            'Isotropic_Phase_Function',
            'Legendre_Phase_Function',
            'Rayleigh_Phase_Function',
            'Reynolds_McCormick_Phase_Function',
            'Table_Phase_Function',
        },
        'SurfaceParticleSizeDistribution_Ptr': {
            'CC1246E_SurfaceParticleSizeDistribution',
            'Regular_SurfaceParticleSizeDistribution',
        },
        'VolumeParticleSizeDistribution_Ptr': {
            'Regular_VolumeParticleSizeDistribution',
        },
        'Distribution_Ptr': {
            'Bimodal_Distribution',
            'Log_Normal_Distribution',
            'Modified_Gamma_Distribution',
            'Normal_Distribution',
            'Weibull_Distribution',
        },
    }

    @staticmethod
    def _set_grating_boundary(node, name, value):
        node["params"][name] = value
        children = node.get("children", {})
        child_name = "grating" if node["model"] == "OneD_CrossGrating" else (
            "top" if name == "medium_i" else "bottom")
        if child_name in children:
            BRDFForm._set_grating_boundary(children[child_name], name, value)

    def _grating_expression_node(self, model, params, *, cross=False, depth=0):
        if depth > 12:
            raise ValueError("Nested grating expressions are too deep.")
        specs = CROSS_GRATING_SPECS if cross else ONE_D_GRATING_SPECS
        if model not in specs:
            raise ValueError(f"Invalid {'cross' if cross else '1D'} grating model: {model}")
        fields = specs[model]["fields"]
        allowed = {field["name"].lower() for field in fields}
        if params.keys() - allowed:
            raise ValueError(f"Unknown parameters for {model}.")
        node = (coerce_cross_grating if cross else coerce_one_d_grating)({"model": model})
        for name in ("medium_i", "medium_t"):
            if name in params:
                self._set_grating_boundary(node, name, params[name])
        for field in fields:
            name = field["name"]
            if name.lower() not in params:
                continue
            value = params[name.lower()]
            if field.get("kind", "param") == "param":
                node["params"][name] = value
            else:
                child_model, child_params = self._parse_model_expression(value)
                for boundary in ("medium_i", "medium_t"):
                    child_params.setdefault(boundary, node["children"][name]["params"][boundary])
                node["children"][name] = self._grating_expression_node(
                    child_model, child_params, cross=field["child_kind"] != "one_d", depth=depth + 1)
        return node

    def _append_model_parameters(self, lines, model, params, *, use_selected_psd=False, depth=0):
        if depth > 12:
            raise ValueError("Nested model expressions are too deep.")
        specs = self.MODEL_PARAM_SPECS.get(model, self.NESTED_MODEL_SPECS.get(model))
        if specs is None:
            raise ValueError(f"Nested model '{model}' is not supported by this input editor.")
        known = {name.lower() for name, *_ in specs}
        supplied = {name.lower(): value for name, value in params.items()}
        unknown = supplied.keys() - known
        if unknown:
            raise ValueError(f"Unknown parameters for {model}: {', '.join(sorted(unknown))}")
        if model in self.MODEL_PARAM_SPECS:
            specs = sorted(specs, key=lambda spec: {"lambda": -3, "substrate": -2, "type": -1}.get(spec[0], 0))
        for name, dtype, default in specs:
            if name == "psd" and use_selected_psd:
                psd_lines = []
                self._append_psd_input_lines(psd_lines)
                for index, val in enumerate(psd_lines[1:], 1):
                    positive = (psd_lines[0] == "Gaussian_PSD_Function" and index == 2) or (
                        psd_lines[0] == "ABC_PSD_Function" and index in (2, 3))
                    finite_float(val, f"{psd_lines[0]} parameter {index}", minimum=0, strict_minimum=positive)
                lines.extend(psd_lines)
                continue
            if model == "CrossRCW_BRDF_Model" and name == "grating" and use_selected_psd:
                tree = self._current_cross_grating_tree()
                errors = validate_cross_grating(tree)
                if errors:
                    raise ValueError("Cross grating validation failed: " + "; ".join(errors))
                substrate = optical_constant(supplied.get("substrate", self.substrate.text()), "Substrate")
                if optical_constant(tree["params"]["medium_t"], "Grating transmitted medium") != substrate:
                    raise ValueError("Cross grating medium_t must match the BRDF substrate.")
                if optical_constant(tree["params"]["medium_i"], "Grating incident medium") != "(1,0)":
                    raise ValueError("CrossRCW BRDF requires grating medium_i to be vacuum (1,0).")
                lines.extend(serialize_cross_grating(tree))
                continue
            value = supplied.get(name.lower(), default)
            if model in {"RCW_BRDF_Model", "CrossRCW_BRDF_Model"}:
                if name == "alpha":
                    finite_float(value, "Detector half angle", minimum=0, strict_minimum=True)
                elif name in {"order", "order1", "order2"}:
                    integer(value, f"Model parameter {name}", minimum=0)
            if dtype in {"double", "int", "dielectric_function", "dielectric_stack", "Table", "string"}:
                lines.append(self._validated_parameter(name, dtype, value))
                continue
            nested_model, nested_params = self._parse_model_expression(single_line(value, f"Model parameter {name}"))
            allowed = self.NESTED_MODEL_TYPES.get(dtype)
            if dtype == "BRDF_Model_Ptr":
                allowed = self.MODEL_PARAM_SPECS
            elif dtype == "Grating_Ptr":
                allowed = ONE_D_GRATING_SPECS
            elif dtype == "CrossGrating_Ptr":
                allowed = CROSS_GRATING_SPECS
            if allowed is None or nested_model not in allowed:
                raise ValueError(f"Model parameter {name} requires {dtype}; {nested_model} is incompatible.")
            if dtype == "BRDF_Model_Ptr":
                defaults = {key.lower(): default_value for key, _, default_value in specs}
                for shared in ("lambda", "substrate", "type"):
                    nested_params.setdefault(shared, supplied.get(shared, defaults[shared]))
            if dtype in {"Free_Space_Scatterer_Ptr", "SphericalScatterer_Ptr"}:
                nested_params.setdefault("lambda", supplied.get("lambda", self.wavelength.text()))
                if model == "Subsurface_Particle_BRDF_Model":
                    nested_params.setdefault("medium", supplied.get("substrate", self.substrate.text()))
            if dtype in {"Grating_Ptr", "CrossGrating_Ptr"}:
                cross = dtype == "CrossGrating_Ptr"
                if model in {"RCW_BRDF_Model", "CrossRCW_BRDF_Model"}:
                    substrate = optical_constant(supplied.get("substrate", self.substrate.text()), "Substrate")
                    nested_params.setdefault("medium_t", substrate)
                    if optical_constant(nested_params["medium_t"], "Grating transmitted medium") != substrate:
                        raise ValueError("Grating medium_t must match the BRDF substrate.")
                node = self._grating_expression_node(nested_model, nested_params, cross=cross, depth=depth)
                if model == "CrossRCW_BRDF_Model" and optical_constant(
                        node["params"]["medium_i"], "Grating incident medium") != "(1,0)":
                    raise ValueError("CrossRCW BRDF requires grating medium_i to be vacuum (1,0).")
                errors = (validate_cross_grating if cross else validate_one_d_grating)(node)
                if errors:
                    raise ValueError("Grating validation failed: " + "; ".join(errors))
                lines.extend((serialize_cross_grating if cross else serialize_one_d_grating)(node))
            else:
                lines.append(nested_model)
                self._append_model_parameters(lines, nested_model, nested_params, depth=depth + 1)
        if model == "Polydisperse_Sphere_BRDF_Model":
            start = finite_float(supplied.get("dstart", "0.1"), "Starting diameter")
            end = finite_float(supplied.get("dend", "100"), "Ending diameter")
            if end < start:
                raise ValueError("Ending diameter must be at least starting diameter.")
            step = finite_float(supplied.get("dstep", "0.01"), "Fractional diameter step")
            if any(bound * (1.0 + step) <= bound for bound in (start, end)):
                raise ValueError("Fractional diameter step is too small to advance the size distribution.")


    def run_brdfprog(self):
        if getattr(self, "_solver_running", False):
            self.output_box.append("A BRDF calculation is already running.")
            return
        self.last_input_path = self.last_stdout_path = self.last_csv_path = self.last_output_meta = None
        self.clear_plot()
        try:
            input_lines = self._build_input_lines()
            exe = shutil.which("brdfprog")
            if not exe:
                raise ValueError(f"'brdfprog' not found. Install SCATMECH or configure SCATMECH_BIN ({scatmech_bin}).")
            timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S_%f")
            data_dir = get_data_dir()
            data_dir.mkdir(parents=True, exist_ok=True)
            input_filename = str(data_dir / f"brdf_input_{timestamp}.txt")
            output_filename = str(data_dir / f"brdf_output_{timestamp}.txt")
            csv_filename = str(data_dir / f"brdf_output_{timestamp}.csv")
            deck = "\n".join(input_lines) + "\n"
            Path(input_filename).write_text(deck, encoding="utf-8")
            self.last_input_path = input_filename
            self.output_box.setPlainText(f"Saved input deck: {input_filename}\nRunning BRDFProg: {exe}")
            meta = self._build_output_meta(csv_filename=csv_filename, stdout_path=output_filename,
                                          input_path=input_filename, timestamp=timestamp,
                                          header_tokens=self.OUTPUT_COLUMNS, column_count=len(self.OUTPUT_COLUMNS))
            result = run_solver(exe, input=deck, parent=self, cwd=data_dir)
            raw_log = result.stdout
            if result.stderr:
                raw_log += "\n# STDERR\n" + "\n".join("# " + line for line in result.stderr.splitlines()) + "\n"
            Path(output_filename).write_text(raw_log, encoding="utf-8")
            self.last_stdout_path = output_filename
            if result.returncode != 0 or "SCATMECH: Error" in result.stderr:
                detail = result.stderr[result.stderr.find("SCATMECH: Error"):] if "SCATMECH: Error" in result.stderr else result.stderr
                raise ValueError(f"brdfprog failed (exit {result.returncode}):\n{detail}")
            try:
                numeric_rows = parse_brdf_output(result.stdout)
            except ValueError as exc:
                raise ValueError(f"{exc}\n{result.stderr.strip()}") from exc
            with open(csv_filename, "w", newline="", encoding="utf-8") as csvfile:
                csv.writer(csvfile).writerows(numeric_rows)
            Path(csv_filename + ".meta.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
            self.last_csv_path = csv_filename
            self.last_output_meta = meta
            self.output_box.append(f"Saved CSV: {csv_filename}")
            self.render_with_external(csv_filename)
        except SolverCancelled:
            self.output_box.append("BRDF run cancelled.")
        except (ValueError, OSError, RuntimeError) as exc:
            self.output_box.setPlainText(f"[Error] {exc}")

    def _safe_float(self, text):
        if text is None:
            return None
        try:
            return float(str(text).strip())
        except Exception:
            return None

    def _build_output_meta(
        self,
        *,
        csv_filename: str,
        stdout_path: str,
        input_path: str,
        timestamp: str,
        header_tokens,
        column_count: int,
    ) -> dict:
        scatter = {
            "incident": self._safe_float(self.incident_angle.text()),
            "scatter_start": self._safe_float(self.scatter_start.text()),
            "scatter_end": self._safe_float(self.scatter_end.text()),
            "scatter_step": self._safe_float(self.scatter_step.text()),
            "azimuth_start": self._safe_float(self.azimuth_start.text()),
            "azimuth_end": self._safe_float(self.azimuth_end.text()),
            "azimuth_step": self._safe_float(self.azimuth_step.text()),
        }

        meta = {
            "csv_path": csv_filename,
            "stdout_path": stdout_path,
            "input_path": input_path,
            "timestamp": timestamp,
            "column_count": column_count,
            "header_tokens": [str(tok) for tok in header_tokens] if header_tokens else [],
            "scatter": scatter,
            "model_family": self.family_selector.currentText(),
            "model_name": self.subclass_selector.currentText() or self.family_selector.currentText(),
            "psd_function": self.psd_function.currentText(),
            "direction_label": self.direction.currentText(),
            "incident_polarization": {"label": "p", "stokes": [1, -1, 0, 0]},
        }

        wavelength_val = self._safe_float(self.wavelength.text())
        if wavelength_val is not None:
            meta["wavelength_um"] = wavelength_val
        substrate_text = self.substrate.text().strip() if hasattr(self, "substrate") else ""
        if substrate_text:
            meta["substrate"] = substrate_text
        meta["direction_code"] = self.DIRECTION_CODES.get(self.direction.currentText(), "0")
        return meta

    def _load_output_meta(self, csv_path: str):
        if not csv_path:
            return None
        try:
            meta_path = Path(str(csv_path) + ".meta.json")
            if not meta_path.exists():
                return None
            meta = json.loads(meta_path.read_text(encoding="utf-8"))
            return meta if isinstance(meta, dict) else None
        except Exception:
            return None


    def _get_text(self, name: str, default: str = "") -> str:
        w = getattr(self, name, None)
        return w.text().strip() if w is not None else default

    def _set_text(self, name: str, value):
        w = getattr(self, name, None)
        if w is not None and value is not None:
            w.setText(str(value))

    def to_params(self) -> dict:
        self._ensure_custom_profiles()
        documents = {}
        if self._current_model_name() == "CrossRCW_BRDF_Model" and self.cross_grating_editor is not None:
            documents = editor_profiles(self.cross_grating_editor)
        elif self._current_model_name() == "RCW_BRDF_Model" and self.custom_profile_binding.document is not None:
            documents = {"grating": self.custom_profile_binding.saved()}
        data = {
            "custom_profiles": documents,
            "program": "BRDFProg",
            "angles": {
                "incident_deg": self._get_text("incident_angle"),
                "scatter_start_deg": self._get_text("scatter_start"),
                "scatter_end_deg": self._get_text("scatter_end"),
                "scatter_step_deg": self._get_text("scatter_step"),
                "azimuth_start_deg": self._get_text("azimuth_start"),
                "azimuth_end_deg": self._get_text("azimuth_end"),
                "azimuth_step_deg": self._get_text("azimuth_step"),
            },
            "general": {
                "wavelength_um": self._get_text("wavelength"),
                "substrate": self._get_text("substrate"),
                "direction": self.direction.currentText() if hasattr(self, "direction") else "",
            },
            "family": self.family_selector.currentText() if hasattr(self, "family_selector") else "",
            "model": self.subclass_selector.currentText() if hasattr(self, "subclass_selector") else "",
            "psd_active": self.psd_function.currentText() if hasattr(self, "psd_function") else "",
            "psd_all": {},
        }
        data["model_params"] = self._collect_model_params()
        if data["model"] == "CrossRCW_BRDF_Model":
            data["grating_tree"] = self._current_cross_grating_tree()


        data["psd_all"]["ABC_PSD_Function"] = {
            "A_um4": self._get_text("psd_A"),
            "B_um": self._get_text("psd_B"),
            "C": self._get_text("psd_C"),
        }
        data["psd_all"]["Fractal_PSD_Function"] = {
            "A_um4": self._get_text("psd_fractalAmp"),
            "gamma": self._get_text("psd_fractalExp"),
        }
        data["psd_all"]["Gaussian_PSD_Function"] = {
            "sigma_um": self._get_text("psd_sigma"),
            "Lc_um": self._get_text("psd_lc"),
        }
        data["psd_all"]["Elliptical_Mesa_PSD_Function"] = {
            "axisx_um": self._get_text("psd_ellipticalX"),
            "axisy_um": self._get_text("psd_ellipticalY"),
            "height_um": self._get_text("psd_mesaHeight"),
            "density_um^-2": self._get_text("psd_mesaDensity"),
        }
        data["psd_all"]["Rectangular_Mesa_PSD_Function"] = {
            "lengthx_um": self._get_text("psd_rectLenX"),
            "lengthy_um": self._get_text("psd_rectLenY"),
            "height_um": self._get_text("psd_rectHeight"),
            "density_um^-2": self._get_text("psd_rectDensity"),
        }
        data["psd_all"]["Triangular_Mesa_PSD_Function"] = {
            "side_um": self._get_text("psd_triSide"),
            "height_um": self._get_text("psd_triHeight"),
            "density_um^-2": self._get_text("psd_triDensity"),
        }
        data["psd_all"]["Rectangular_Pyramid_PSD_Function"] = {
            "lengthx_um": self._get_text("psd_pyrLenX"),
            "lengthy_um": self._get_text("psd_pyrLenY"),
            "height_um": self._get_text("psd_pyrHeight"),
            "density_um^-2": self._get_text("psd_pyrDensity"),
        }
        data["psd_all"]["Triangular_Pyramid_PSD_Function"] = {
            "side_um": self._get_text("psd_tpSide"),
            "height_um": self._get_text("psd_tpHeight"),
            "density_um^-2": self._get_text("psd_tpDensity"),
        }
        data["psd_all"]["Parabolic_Dimple_PSD_Function"] = {
            "axisx_um": self._get_text("psd_pdAxisX"),
            "axisy_um": self._get_text("psd_pdAxisY"),
            "height_um": self._get_text("psd_pdHeight"),
            "density_um^-2": self._get_text("psd_pdDensity"),
        }
        return data

    def from_params(self, p: dict):
        if not isinstance(p, dict):
            raise ValueError("BRDF parameters must be a JSON object.")
        if not isinstance(p.get("custom_profiles", {}), dict):
            raise ValueError("custom_profiles must be an object.")
        if p.get("program", "BRDFProg") != "BRDFProg":
            raise ValueError("The parameter file is not for BRDFProg.")
        for section in ("angles", "general", "model_params", "psd_all"):
            if not isinstance(p.get(section, {}), dict):
                raise ValueError(f"BRDF {section} must be an object.")
        for section in ("angles", "general", "model_params"):
            if any(isinstance(value, (dict, list)) for value in p.get(section, {}).values()):
                raise ValueError(f"BRDF {section} values must be scalar values.")
        for name, values in p.get("psd_all", {}).items():
            if not isinstance(values, dict) or any(isinstance(value, (dict, list)) for value in values.values()):
                raise ValueError(f"PSD parameters for {name} must be an object of scalar values.")
        family = p.get("family", self.family_selector.currentText())
        if not isinstance(family, str) or self.family_selector.findText(family) < 0:
            raise ValueError(f"Unknown BRDF model family: {family}")
        model = p.get("model")
        if model is not None:
            if not isinstance(model, str) or model not in self.MODEL_PARAM_SPECS:
                raise ValueError(f"Unknown BRDF model: {model}")

            if model not in self.MODEL_FAMILIES.get(family, {family}):
                raise ValueError(f"BRDF model {model} does not belong to {family}.")
        selected_model = model or (self.MODEL_FAMILIES.get(family) or [family])[0]
        if family == self.family_selector.currentText() and model is None:
            selected_model = self._current_model_name()
        model_params = dict(p.get("model_params", {}))
        if "Axisymmetric_Particle_BRDF_Model" in selected_model and "nmax" in model_params:
            if "mmax" in model_params and str(model_params["mmax"]) != str(model_params["nmax"]):
                raise ValueError("Conflicting saved mmax and legacy nmax values.")
            model_params["mmax"] = model_params.pop("nmax")
        allowed = {name for name, *_ in self.MODEL_PARAM_SPECS.get(selected_model, [])}
        unknown = set(model_params) - allowed
        if unknown:
            raise ValueError(f"Unknown parameters for {selected_model}: {', '.join(sorted(unknown))}")
        for name in p.get("psd_all", {}):
            if self.psd_function.findText(name) < 0:
                raise ValueError(f"Unknown PSD function: {name}")
        general = dict(p.get("general", {}))
        for name, key in (("lambda", "wavelength_um"), ("substrate", "substrate"), ("type", "direction")):
            if name not in model_params:
                continue
            value = model_params.pop(name)
            if name == "type":
                code = str(integer(value, "BRDF direction", minimum=0, maximum=3))
                value = next(label for label, candidate in self.DIRECTION_CODES.items() if candidate == code)
            if key in general:
                normalize = (lambda v: optical_constant(v, "Substrate")) if name == "substrate" else (
                    (lambda v: finite_float(v, "Wavelength")) if name == "lambda" else str)
                if normalize(general[key]) != normalize(value):
                    raise ValueError(f"Conflicting saved values for {name}.")
            general[key] = value
        if "direction" in general and general["direction"] not in self.DIRECTION_CODES:
            raise ValueError(f"Unknown BRDF direction: {general['direction']}")
        psd = p.get("psd_active")
        if "psd" in model_params:
            legacy_psd, psd_params = self._parse_model_expression(model_params.pop("psd"))
            if psd_params:
                raise ValueError("Save PSD parameters in psd_all and select the function with psd_active.")
            if psd is not None and psd != legacy_psd:
                raise ValueError("Conflicting saved PSD selections.")
            psd = legacy_psd
        if psd is not None and (not isinstance(psd, str) or self.psd_function.findText(psd) < 0):
            raise ValueError(f"Unknown PSD function: {psd}")
        tree = p.get("grating_tree")
        if selected_model == "CrossRCW_BRDF_Model" and "grating" in model_params:
            if tree is not None:
                raise ValueError("Specify only one cross grating, using grating_tree or model_params.grating.")
            grating_model, grating_params = self._parse_model_expression(model_params.pop("grating"))
            tree = self._grating_expression_node(grating_model, grating_params, cross=True)
        profile_documents = {}
        if p.get("custom_profiles"):
            if selected_model == "CrossRCW_BRDF_Model":
                if tree is None:
                    raise ValueError("Custom profiles require a saved grating_tree.")
                profile_documents, tree = prepare_profiles(p["custom_profiles"], tree)
            elif selected_model == "RCW_BRDF_Model":
                name, parameters = self._parse_model_expression(model_params.get("grating", ""))
                parameters.setdefault("medium_t", general.get("substrate", self.substrate.text()))
                one_d_tree = self._grating_expression_node(name, parameters)
                profile_documents, one_d_tree = prepare_profiles(p["custom_profiles"], one_d_tree)
                model_params["grating"] = node_expression(one_d_tree)
            else:
                raise ValueError("Custom profiles are only supported by RCW and CrossRCW models.")
        if tree is not None:
            errors = validate_cross_grating(tree)
            if errors:
                raise ValueError("Invalid saved grating: " + "; ".join(errors))

        self.custom_profile_binding.restore(None, None)
        self._cross_profile_state = {}
        for key, attr in {"incident_deg": "incident_angle", "scatter_start_deg": "scatter_start",
                          "scatter_end_deg": "scatter_end", "scatter_step_deg": "scatter_step",
                          "azimuth_start_deg": "azimuth_start", "azimuth_end_deg": "azimuth_end",
                          "azimuth_step_deg": "azimuth_step"}.items():
            self._set_text(attr, p.get("angles", {}).get(key))
        self._set_text("wavelength", general.get("wavelength_um"))
        self._set_text("substrate", general.get("substrate"))
        if "direction" in general:
            self.direction.setCurrentText(general["direction"])
        self.family_selector.setCurrentText(family)
        if model is not None:
            self.subclass_selector.setCurrentText(model)
        if psd is not None:
            self.psd_function.setCurrentText(psd)
        psd_fields = {
            "ABC_PSD_Function": {"A_um4": "psd_A", "B_um": "psd_B", "C": "psd_C"},
            "Fractal_PSD_Function": {"A_um4": "psd_fractalAmp", "gamma": "psd_fractalExp"},
            "Gaussian_PSD_Function": {"sigma_um": "psd_sigma", "Lc_um": "psd_lc"},
            "Elliptical_Mesa_PSD_Function": {"axisx_um": "psd_ellipticalX", "axisy_um": "psd_ellipticalY",
                "height_um": "psd_mesaHeight", "density_um^-2": "psd_mesaDensity"},
            "Rectangular_Mesa_PSD_Function": {"lengthx_um": "psd_rectLenX", "lengthy_um": "psd_rectLenY",
                "height_um": "psd_rectHeight", "density_um^-2": "psd_rectDensity"},
            "Triangular_Mesa_PSD_Function": {"side_um": "psd_triSide", "height_um": "psd_triHeight", "density_um^-2": "psd_triDensity"},
            "Rectangular_Pyramid_PSD_Function": {"lengthx_um": "psd_pyrLenX", "lengthy_um": "psd_pyrLenY",
                "height_um": "psd_pyrHeight", "density_um^-2": "psd_pyrDensity"},
            "Triangular_Pyramid_PSD_Function": {"side_um": "psd_tpSide", "height_um": "psd_tpHeight", "density_um^-2": "psd_tpDensity"},
            "Parabolic_Dimple_PSD_Function": {"axisx_um": "psd_pdAxisX", "axisy_um": "psd_pdAxisY",
                "height_um": "psd_pdHeight", "density_um^-2": "psd_pdDensity"},
        }
        for name, values in p.get("psd_all", {}).items():
            for key, attr in psd_fields.get(name, {}).items():
                self._set_text(attr, values.get(key))
        for name, value in model_params.items():
            if name in self.param_widgets and value is not None:
                self.param_widgets[name].setText(str(value))
        if tree is not None and self.cross_grating_editor is not None:
            self.cross_grating_editor.from_node(tree)
            self._cross_grating_state = self.cross_grating_editor.to_node()
            restore_editor_profiles(self.cross_grating_editor, profile_documents)
            self._cross_profile_state = profile_documents
        elif selected_model == "RCW_BRDF_Model":
            self.custom_profile_binding.restore(profile_documents.get("grating"), self._custom_profile_node())
        return True


    def save_to_json(self, path: str):
        import json, pathlib
        data = self.to_params()
        data.setdefault("_schema", {"program": "BRDFProg", "version": 2})
        pathlib.Path(path).write_text(json.dumps(data, indent=2), encoding="utf-8")

    def load_from_json(self, path: str):
        import json, pathlib
        data = json.loads(pathlib.Path(path).read_text(encoding="utf-8"))
        self.from_params(data)


    def _on_model_changed(self, *_):

        try:
            self.populate_model_params()
        except Exception as e:
            print("Param populate error:", e)

    def open_last_output(self):

        from brdfplot import _load_numeric_table
        from PyQt5.QtWidgets import QDialog, QTableWidget, QTableWidgetItem

        data_dir = get_data_dir()
        candidates = [getattr(self, "last_csv_path", None), getattr(self, "last_stdout_path", None)]
        if not any(path and Path(path).is_file() for path in candidates):
            candidates = sorted(data_dir.glob("brdf_output_*"), key=lambda path: path.stat().st_mtime, reverse=True)
        for path in candidates:
            if not path or Path(path).suffix not in {".csv", ".txt"} or not Path(path).is_file():
                continue
            try:
                data = _load_numeric_table(str(path))
            except (OSError, ValueError) as exc:
                self.output_box.append(f"Could not read {path}: {exc}")
                continue
            ncols = data.shape[1]
            header = self.OUTPUT_COLUMNS if ncols == 7 else [f"C{i + 1}" for i in range(ncols)]
            dlg = QDialog(self)
            dlg.setWindowTitle(f"Output Table: {Path(path).name}")
            layout = QVBoxLayout(dlg)
            table = QTableWidget(data.shape[0], ncols, dlg)
            table.setHorizontalHeaderLabels(header)
            for row, values in enumerate(data):
                for column, value in enumerate(values):
                    table.setItem(row, column, QTableWidgetItem(f"{value:g}"))
            layout.addWidget(table)
            dlg.resize(900, 600)
            dlg.exec_()
            return
        self.output_box.append(f"No readable BRDF output file found in {data_dir}.")

    def open_last_input(self):

        data_dir = str(get_data_dir())
        path = getattr(self, "last_input_path", None)
        if not path or not os.path.exists(path):
            try:
                cand = [os.path.join(data_dir, f) for f in os.listdir(data_dir) if f.startswith("brdf_input_") and f.endswith(".txt")]
                path = max(cand, key=os.path.getmtime) if cand else None
            except Exception:
                path = None
        if not path:
            self.output_box.append(f"No BRDF input file found in {data_dir}.")
            return
        from PyQt5.QtWidgets import QDialog, QVBoxLayout, QTextEdit
        dlg = QDialog(self)
        dlg.setWindowTitle(f"Input: {os.path.basename(path)}")
        layout = QVBoxLayout(dlg)
        txt = QTextEdit(dlg)
        txt.setReadOnly(True)
        try:
            txt.setPlainText(Path(path).read_text(encoding="utf-8"))
        except OSError as exc:
            self.output_box.append(f"Could not open input: {exc}")
            return
        layout.addWidget(txt)
        dlg.resize(700, 500)
        dlg.exec_()
