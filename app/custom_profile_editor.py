from copy import deepcopy
from dataclasses import asdict
import json
import math
import uuid

import numpy as np
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QKeySequence
from PyQt5.QtWidgets import (
    QApplication, QComboBox, QDoubleSpinBox, QFileDialog, QFormLayout,
    QGroupBox, QHBoxLayout, QLabel, QLineEdit, QMessageBox, QPushButton,
    QScrollArea, QShortcut, QSizePolicy, QSpinBox, QSplitter, QStackedLayout,
    QStackedWidget, QTabWidget, QVBoxLayout, QWidget,
)

from custom_profile import default_editor_profile, default_profile
from custom_profile_canvas import move_profile_point
from custom_profile_constraints import (
    enforce_region_constraints, ensure_spacer_constraints, region_delete_error,
    region_vertex_edit_error, remap_region_vertices,
)
from custom_profile_regions_canvas import CustomProfileRegionsCanvas, move_region_vertex, region_mode
from custom_profile_1d_controls import NativeDefinitionPanel, RegionProperties, region_base_layer
from custom_profile_preview import CustomProfilePreview


def document_signature(document):
    return json.dumps(asdict(document), sort_keys=True, separators=(",", ":"))


class CustomProfileEditor(QWidget):
    def __init__(self, parent=None, document=None, context_provider=None, apply_callback=None):
        super().__init__(parent, Qt.Window)
        self.setWindowModality(Qt.WindowModal)
        self.setAttribute(Qt.WA_DeleteOnClose, False)
        self.setWindowTitle("Custom Grating Profile")
        self.setMinimumSize(900, 600)
        screen = self.screen() or QApplication.primaryScreen()
        available = screen.availableGeometry() if screen is not None else None
        self.resize(min(1560, max(900, available.width() - 48)) if available else 1560,
                    min(900, max(600, available.height() - 80)) if available else 900)
        self.context_provider = context_provider or (lambda: {})
        self.apply_callback = apply_callback
        self.document = deepcopy(document if document is not None else default_editor_profile(1))
        ensure_spacer_constraints(self.document)
        enforce_region_constraints(self.document)
        self._applied_document = deepcopy(self.document)
        self._applied_signature = document_signature(self.document)
        self._history = [deepcopy(self.document)]
        self._history_index = 0
        self._syncing = False
        self._drag_before = None
        self._shown = False
        self._fields = {}
        self._field_rows = {}
        self._build_ui()
        self._sync_controls()
        self._update_preview()
        application = QApplication.instance()
        if application is not None:
            application.aboutToQuit.connect(self.shutdown)

    @property
    def dirty(self):
        return document_signature(self.document) != self._applied_signature

    def _build_ui(self):
        layout = QVBoxLayout(self)
        header = QHBoxLayout()
        self.undo_button = QPushButton("Undo")
        self.undo_button.clicked.connect(self.undo)
        self.redo_button = QPushButton("Redo")
        self.redo_button.clicked.connect(self.redo)
        self.reset_button = QPushButton("Reset to default")
        self.reset_button.clicked.connect(self.reset)
        self.fit_button = QPushButton("Fit drawing")
        self.fit_button.clicked.connect(lambda: self.canvas.fit())
        for button in (self.undo_button, self.redo_button, self.reset_button, self.fit_button):
            header.addWidget(button)
        header.addStretch()
        self.draft_label = QLabel()
        header.addWidget(self.draft_label)
        layout.addLayout(header)
        self.instructions = QLabel()
        self.instructions.setWordWrap(True)
        layout.addWidget(self.instructions)
        self.one_d_controls = QWidget()
        mode_row = QHBoxLayout(self.one_d_controls)
        mode_row.setContentsMargins(0, 0, 0, 0)
        mode_row.addWidget(QLabel("1D editing mode"))
        self.one_d_mode = QComboBox()
        for label, value in (("Surface profile", "surface"), ("Material regions", "regions"), ("Native definition", "native")):
            self.one_d_mode.addItem(label, value)
        self.one_d_mode.currentIndexChanged.connect(lambda: self.set_one_d_mode(self.one_d_mode.currentData()))
        mode_row.addWidget(self.one_d_mode)
        self.spacer_preset_button = QPushButton("Load core + spacers example")
        self.spacer_preset_button.clicked.connect(self.load_spacer_preset)
        mode_row.addWidget(self.spacer_preset_button)
        mode_row.addStretch()
        layout.addWidget(self.one_d_controls)
        self.splitter = QSplitter(Qt.Horizontal)
        self.visual_container = QWidget()
        self.visual_stack = QStackedLayout(self.visual_container)
        self.visual_splitter = QSplitter(Qt.Horizontal)
        self.visual_tabs = QTabWidget()
        self.visual_stack.addWidget(self.visual_splitter)
        self.visual_stack.addWidget(self.visual_tabs)
        self.splitter.addWidget(self.visual_container)
        self._compact_layout = None
        drawing = QWidget()
        self._drawing_panel = drawing
        drawing_layout = QVBoxLayout(drawing)
        drawing_layout.setContentsMargins(0, 0, 0, 0)
        self.canvas = CustomProfileRegionsCanvas(self)
        self.canvas.selection_changed.connect(self._selection_changed)
        self.canvas.edit_started.connect(self._begin_drag)
        self.canvas.geometry_edited.connect(self._drag_changed)
        self.canvas.edit_finished.connect(self._finish_drag)
        self.canvas.polygon_created.connect(self.add_drawn_region)
        self.drawing_stack = QStackedWidget()
        self.drawing_stack.addWidget(self.canvas)
        self.native_panel = NativeDefinitionPanel(self)
        self.drawing_stack.addWidget(self.native_panel)
        drawing_layout.addWidget(self.drawing_stack)
        self.drawing_actions = QHBoxLayout()
        self.add_point_button = QPushButton("Add point")
        self.add_point_button.clicked.connect(self.add_point)
        self.delete_point_button = QPushButton("Delete point")
        self.delete_point_button.clicked.connect(self.delete_point)
        self.shape_selector = QComboBox()
        for label, key in (("Circle", "circle"), ("Rectangle", "rectangle"), ("Polygon", "polygon")):
            self.shape_selector.addItem(label, key)
        self.add_feature_button = QPushButton("Add feature")
        self.add_feature_button.clicked.connect(lambda: self.add_feature(self.shape_selector.currentData()))
        self.duplicate_button = QPushButton("Duplicate")
        self.duplicate_button.clicked.connect(self.duplicate_feature)
        self.delete_feature_button = QPushButton("Delete")
        self.delete_feature_button.clicked.connect(self.delete_feature)
        for widget in (self.add_point_button, self.delete_point_button, self.shape_selector, self.add_feature_button, self.duplicate_button, self.delete_feature_button):
            self.drawing_actions.addWidget(widget)
        drawing_layout.addLayout(self.drawing_actions)
        self.region_actions = QWidget()
        region_actions = QHBoxLayout(self.region_actions)
        region_actions.setContentsMargins(0, 0, 0, 0)
        self.region_shape = QComboBox()
        for label, value in (("Rectangle", "rectangle"), ("Trapezoid", "trapezoid"), ("Ellipse", "ellipse"), ("Full-period film", "film")):
            self.region_shape.addItem(label, value)
        self.add_region_button = QPushButton("Add region")
        self.add_region_button.clicked.connect(lambda: self.add_region(self.region_shape.currentData()))
        self.draw_region_button = QPushButton("Draw polygon")
        self.draw_region_button.clicked.connect(self.begin_region_drawing)
        self.duplicate_region_button = QPushButton("Duplicate")
        self.duplicate_region_button.clicked.connect(self.duplicate_region)
        self.delete_region_button = QPushButton("Delete")
        self.delete_region_button.clicked.connect(self.delete_region)
        for widget in (self.region_shape, self.add_region_button, self.draw_region_button, self.duplicate_region_button, self.delete_region_button):
            region_actions.addWidget(widget)
        drawing_layout.addWidget(self.region_actions)
        self.visual_splitter.addWidget(drawing)
        self.preview = CustomProfilePreview(parent=self)
        self.preview.validation_changed.connect(self._validation_changed)
        self.visual_splitter.addWidget(self.preview)
        scroll = QScrollArea()
        self.property_scroll = scroll
        scroll.setWidgetResizable(True)
        scroll.setMinimumWidth(290)
        scroll.setMaximumWidth(390)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        properties = QWidget()
        self.properties_layout = QVBoxLayout(properties)
        geometry_group = QGroupBox("Geometry (µm)")
        geometry_form = QFormLayout(geometry_group)
        geometry_form.setRowWrapPolicy(QFormLayout.WrapLongRows)
        for key, label, lo, hi in (
            ("period_um", "Period", .000001, 1000000),
            ("d1_um", "Period a₁", .000001, 1000000),
            ("d2_um", "Period a₂", .000001, 1000000),
            ("zeta_deg", "Skew (degrees)", -89.99, 89.99),
            ("film_thickness_um", "Film thickness", .000001, 1000000),
        ):
            field = self._number(geometry_form, key, label, lo, hi)
            field.valueChanged.connect(lambda value, name=key: self.set_parameter(name, value))
        self.film_material = self._material(geometry_form, "Film material", lambda value: self.set_parameter("film_material", value))
        self.film_material_label = geometry_form.labelForField(self.film_material.parentWidget())
        self.media_label = QLabel()
        self.media_label.setWordWrap(True)
        geometry_form.addRow("Inherited media", self.media_label)
        self.properties_layout.addWidget(geometry_group)
        self.selection_group = QGroupBox("Selection")
        selection_form = QFormLayout(self.selection_group)
        selection_form.setRowWrapPolicy(QFormLayout.WrapLongRows)
        self.selection_combo = QComboBox()
        self.selection_combo.setMinimumWidth(100)
        self.selection_combo.setSizePolicy(QSizePolicy.Ignored, QSizePolicy.Fixed)
        self.selection_combo.currentIndexChanged.connect(lambda index: self._choose_selection(index))
        selection_form.addRow("Selected", self.selection_combo)
        self.feature_name = QLineEdit()
        self.feature_name.setMinimumWidth(50)
        self.feature_name.setSizePolicy(QSizePolicy.Ignored, QSizePolicy.Fixed)
        self.feature_name.editingFinished.connect(lambda: self.set_feature_parameter("name", self.feature_name.text()))
        self._add_row(selection_form, "name", "Name", self.feature_name)
        self.role_combo = QComboBox()
        self.role_combo.addItem("Hole", "hole")
        self.role_combo.addItem("Post", "post")
        self.role_combo.currentIndexChanged.connect(lambda: self.set_feature_parameter("role", self.role_combo.currentData()))
        self._add_row(selection_form, "role", "Feature type", self.role_combo)
        for key, label in (("point_x", "Point x"), ("point_z", "Point z"), ("x_um", "Center x"), ("y_um", "Center y"), ("radius_um", "Radius"), ("width_um", "Width"), ("height_um", "Height in y"), ("rotation_deg", "Rotation (degrees)"), ("depth_um", "Hole depth"), ("post_height_um", "Post height")):
            field = self._number(selection_form, key, label)
            if key.startswith("point_"):
                field.valueChanged.connect(lambda value, name=key: self.set_point_coordinate(name, value))
            else:
                field.valueChanged.connect(lambda value, name=key: self.set_feature_parameter(name, value))
        self.post_material = self._material(selection_form, "Post material", lambda value: self.set_feature_parameter("material", value), "post_material")
        self.vertex_combo = QComboBox()
        self.vertex_combo.currentIndexChanged.connect(self._choose_vertex)
        self._add_row(selection_form, "vertex_combo", "Polygon vertex", self.vertex_combo)
        for key, label in (("vertex_x", "Vertex local x"), ("vertex_y", "Vertex local y")):
            field = self._number(selection_form, key, label)
            field.valueChanged.connect(lambda value, name=key: self.set_vertex_coordinate(name, value))
        vertex_buttons = QWidget()
        vertex_actions = QHBoxLayout(vertex_buttons)
        vertex_actions.setContentsMargins(0, 0, 0, 0)
        self.add_vertex_button = QPushButton("Add vertex")
        self.delete_vertex_button = QPushButton("Delete vertex")
        self.add_vertex_button.clicked.connect(self.add_vertex)
        self.delete_vertex_button.clicked.connect(self.delete_vertex)
        vertex_actions.addWidget(self.add_vertex_button)
        vertex_actions.addWidget(self.delete_vertex_button)
        self._add_row(selection_form, "vertex_buttons", "", vertex_buttons)
        self.properties_layout.addWidget(self.selection_group)
        self.region_properties = RegionProperties(self)
        self.properties_layout.addWidget(self.region_properties)
        solver_group = QGroupBox("Solver and approximation")
        solver_form = QFormLayout(solver_group)
        solver_form.setRowWrapPolicy(QFormLayout.WrapLongRows)
        for key, label, lo, hi in (("nlayers", "Layers", 1, 10000), ("grid1", "Grid a₁", 2, 16384), ("grid2", "Grid a₂", 2, 16384), ("curve_segments", "Circle segments", 12, 2048)):
            field = QSpinBox()
            field.setRange(lo, hi)
            field.setKeyboardTracking(False)
            self._fields[key] = field
            self._add_row(solver_form, key, label, field)
            field.valueChanged.connect(lambda value, name=key: self.set_parameter(name, value))
        self.properties_layout.addWidget(solver_group)
        self.properties_layout.addStretch()
        scroll.setWidget(properties)
        self.splitter.addWidget(scroll)
        self.splitter.setStretchFactor(0, 1)
        self.splitter.setStretchFactor(1, 0)
        self.splitter.setSizes([1100, 300])
        self.visual_splitter.setStretchFactor(0, 4)
        self.visual_splitter.setStretchFactor(1, 5)
        self.visual_splitter.setSizes([510, 650])
        self._adapt_layout()
        layout.addWidget(self.splitter, 1)
        self.status_label = QLabel("Preparing preview…")
        self.status_label.setWordWrap(True)
        self.status_label.setTextFormat(Qt.PlainText)
        layout.addWidget(self.status_label)
        bottom = QHBoxLayout()
        bottom.addWidget(QLabel("Apply updates the simulation inputs. It does not start a simulation."))
        bottom.addStretch()
        self.apply_button = QPushButton("Apply to simulation")
        self.apply_button.setEnabled(False)
        self.apply_button.clicked.connect(self.apply)
        self.close_button = QPushButton("Close")
        self.close_button.clicked.connect(self.close)
        bottom.addWidget(self.apply_button)
        bottom.addWidget(self.close_button)
        layout.addLayout(bottom)
        self._undo_shortcut = QShortcut(QKeySequence.Undo, self, self.undo)
        self._redo_shortcut = QShortcut(QKeySequence.Redo, self, self.redo)

    def _add_row(self, form, key, label, field):
        caption = QLabel(label)
        form.addRow(caption, field)
        self._field_rows[key] = (caption, field)

    def _number(self, form, key, label, lo=-1000000, hi=1000000):
        field = QDoubleSpinBox()
        field.setDecimals(8)
        field.setRange(lo, hi)
        field.setSingleStep(.01)
        field.setKeyboardTracking(False)
        field.setMinimumWidth(80)
        field.setSizePolicy(QSizePolicy.Ignored, QSizePolicy.Fixed)
        self._fields[key] = field
        self._add_row(form, key, label, field)
        return field

    def _material(self, form, label, callback, key=None):
        container = QWidget()
        row = QHBoxLayout(container)
        row.setContentsMargins(0, 0, 0, 0)
        entry = QLineEdit()
        entry.setMinimumWidth(50)
        entry.setSizePolicy(QSizePolicy.Ignored, QSizePolicy.Fixed)
        entry.setPlaceholderText("(n,k) or dielectric file")
        entry.setToolTip("Optical constants (n,k), or a SCATMECH dielectric table file")
        entry.editingFinished.connect(lambda: callback(entry.text()))
        browse = QPushButton("…")
        browse.setFixedWidth(28)
        browse.clicked.connect(lambda: self._browse_material(entry, callback))
        row.addWidget(entry)
        row.addWidget(browse)
        if key:
            self._add_row(form, key, label, container)
        else:
            form.addRow(label, container)
        return entry

    def _browse_material(self, entry, callback):
        path, _ = QFileDialog.getOpenFileName(self, "Select dielectric material", "", "Material tables (*.nk *.dat *.txt);;All files (*)")
        if path:
            entry.setText(path)
            callback(path)

    def _row_visible(self, key, visible):
        for widget in self._field_rows[key]:
            widget.setVisible(visible)

    def _sync_controls(self):
        self._syncing = True
        document = self.document
        one_d = document.dimension == 1
        mode = getattr(document, "one_d_mode", "surface") if one_d else "surface"
        self.one_d_controls.setVisible(one_d)
        self.one_d_mode.setCurrentIndex(max(0, self.one_d_mode.findData(mode)))
        self.region_actions.setVisible(one_d and mode == "regions")
        self.region_properties.setVisible(one_d and mode == "regions")
        self.selection_group.setVisible(not one_d or mode == "surface")
        self.property_scroll.setMinimumWidth(350 if one_d and mode == "regions" else 290)
        self.drawing_stack.setCurrentWidget(self.native_panel if one_d and mode == "native" else self.canvas)
        self.fit_button.setEnabled(not one_d or mode != "native")
        self.film_material.parentWidget().setVisible(not one_d or mode == "surface")
        if self.film_material_label is not None:
            self.film_material_label.setVisible(not one_d or mode == "surface")
        if one_d and mode in ("regions", "native"):
            self._sync_advanced_controls(mode)
            return
        self.instructions.setText("Drag the blue profile points or edit their coordinates. Endpoints repeat periodically; equal neighboring x coordinates make a vertical step. Scroll to zoom." if one_d else "Select and drag a feature. Square handles resize circles/rectangles or move polygon vertices; the diamond rotates. Holes and posts must remain separate and inside the unit cell. Scroll to zoom.")
        self._field_rows["curve_segments"][0].setText("Circle segments")
        for key in ("period_um",):
            self._row_visible(key, one_d)
        for key in ("d1_um", "d2_um", "zeta_deg", "film_thickness_um", "grid1", "grid2", "curve_segments"):
            self._row_visible(key, not one_d)
        for key in ("period_um", "d1_um", "d2_um", "zeta_deg", "film_thickness_um", "nlayers", "grid1", "grid2", "curve_segments"):
            self._fields[key].setValue(getattr(document, key))
        self.film_material.setText(document.film_material)
        self.media_label.setText(f"Incident: {document.medium_i}\nSubstrate: {document.medium_t}")
        self.media_label.setTextFormat(Qt.PlainText)
        self.selection_combo.clear()
        if one_d:
            self.selection_combo.addItems([f"Point {index + 1}" for index in range(len(document.points))])
        else:
            self.selection_combo.addItems([f"{index + 1}. {feature.get('name', feature['shape'].title())}" for index, feature in enumerate(document.features)])
        count = self.selection_combo.count()
        self.canvas.selected_index = min(max(self.canvas.selected_index, 0), count - 1)
        self.selection_combo.setCurrentIndex(self.canvas.selected_index)
        index = self.canvas.selected_index
        selected = index >= 0
        self.selection_group.setEnabled(selected)
        self.add_point_button.setVisible(one_d)
        self.delete_point_button.setVisible(one_d)
        self.delete_point_button.setEnabled(one_d and 0 < index < len(document.points) - 1)
        for widget in (self.shape_selector, self.add_feature_button, self.duplicate_button, self.delete_feature_button):
            widget.setVisible(not one_d)
        self.duplicate_button.setEnabled(selected)
        self.delete_feature_button.setEnabled(selected)
        for key in ("point_x", "point_z"):
            self._row_visible(key, one_d)
        if one_d and selected:
            self._fields["point_x"].setValue(document.points[index][0])
            self._fields["point_z"].setValue(document.points[index][1])
            self._fields["point_x"].setEnabled(index not in (0, len(document.points) - 1))
        feature = document.features[index] if not one_d and selected else None
        shape = feature["shape"] if feature else None
        role = feature["role"] if feature else None
        visible = {"name": not one_d, "role": not one_d, "x_um": not one_d, "y_um": not one_d, "radius_um": shape == "circle", "width_um": shape == "rectangle", "height_um": shape == "rectangle", "rotation_deg": shape in ("rectangle", "polygon"), "depth_um": role == "hole", "post_height_um": role == "post", "post_material": role == "post"}
        for key, value in visible.items():
            self._row_visible(key, value)
        if feature:
            self.feature_name.setText(feature.get("name", ""))
            self.role_combo.setCurrentIndex(0 if role == "hole" else 1)
            for key in ("x_um", "y_um", "radius_um", "width_um", "height_um", "rotation_deg", "depth_um", "post_height_um"):
                self._fields[key].setValue(feature.get(key, .1 if key not in ("x_um", "y_um", "rotation_deg") else 0))
            self.post_material.setText(feature.get("material", document.film_material))
        for key in ("vertex_combo", "vertex_x", "vertex_y", "vertex_buttons"):
            self._row_visible(key, shape == "polygon")
        self.vertex_combo.clear()
        if shape == "polygon":
            self.vertex_combo.addItems([f"Vertex {i + 1}" for i in range(len(feature["vertices"]))])
            self.canvas.selected_vertex = min(max(self.canvas.selected_vertex, 0), len(feature["vertices"]) - 1)
            self.vertex_combo.setCurrentIndex(self.canvas.selected_vertex)
            vertex = feature["vertices"][self.canvas.selected_vertex]
            self._fields["vertex_x"].setValue(vertex[0])
            self._fields["vertex_y"].setValue(vertex[1])
            self.delete_vertex_button.setEnabled(len(feature["vertices"]) > 3)
        self.undo_button.setEnabled(self._history_index > 0)
        self.redo_button.setEnabled(self._history_index < len(self._history) - 1)
        self.draft_label.setText("Unapplied changes" if self.dirty else "No unapplied edits")
        self._syncing = False
        self.canvas.set_document(document)

    def _sync_advanced_controls(self, mode):
        document = self.document
        self.instructions.setText("Drag a region or its vertices; edit exact coordinates and materials on the right. Later regions replace earlier material. Overhangs, shared edges, buried voids, and periodic crossings are supported." if mode == "regions" else "Edit the complete native SCATMECH boundary definition and parameter values. Native definitions are not reconstructed as surface or region shapes when switching modes. Undo preserves access to earlier native edits.")
        self._row_visible("period_um", mode == "regions")
        self._row_visible("nlayers", True)
        for key in ("d1_um", "d2_um", "zeta_deg", "film_thickness_um", "grid1", "grid2"):
            self._row_visible(key, False)
        self._row_visible("curve_segments", mode == "regions")
        self._field_rows["curve_segments"][0].setText("New ellipse segments")
        for key in ("period_um", "nlayers", "curve_segments"):
            self._fields[key].setValue(getattr(document, key))
        self.media_label.setText(f"Incident: {document.medium_i}\nSubstrate: {document.medium_t}")
        self.media_label.setTextFormat(Qt.PlainText)
        for widget in (self.add_point_button, self.delete_point_button, self.shape_selector, self.add_feature_button, self.duplicate_button, self.delete_feature_button):
            widget.hide()
        if mode == "regions":
            self.region_properties.sync()
            selected = self.canvas.selected_index >= 0
            self.duplicate_region_button.setEnabled(selected)
            reason = region_delete_error(document, self.canvas.selected_index) if selected else ""
            self.delete_region_button.setEnabled(selected and not reason)
            self.delete_region_button.setToolTip(reason)
        else:
            self.native_panel.sync()
        self.undo_button.setEnabled(self._history_index > 0)
        self.redo_button.setEnabled(self._history_index < len(self._history) - 1)
        self.draft_label.setText("Unapplied changes" if self.dirty else "No unapplied edits")
        self._syncing = False
        if mode == "regions":
            self.canvas.set_document(document)

    def _choose_selection(self, index):
        if not self._syncing:
            self.canvas.set_selection(index)

    def _choose_vertex(self, index):
        if not self._syncing:
            self.canvas.set_selection(self.canvas.selected_index, index)

    def _selection_changed(self, index, vertex):
        if not self._syncing:
            self._sync_controls()

    def _set_document(self, document, record=True, fit=False):
        enforce_region_constraints(document, previous=self.document)
        if document_signature(document) == document_signature(self.document):
            return
        self.document = deepcopy(document)
        if record:
            self._history = self._history[:self._history_index + 1] + [deepcopy(document)]
            self._history_index += 1
        self._sync_controls()
        if fit:
            self.canvas.fit()
        self._update_preview()

    def _update_preview(self):
        self.apply_button.setEnabled(False)
        try:
            context = self.context_provider()
            self.preview.set_document(deepcopy(self.document), context)
        except Exception as exc:
            self.status_label.setText(str(exc))

    def _validation_changed(self, revision, valid, message):
        if revision != self.preview.revision:
            return
        self.status_label.setText(message)
        self.apply_button.setEnabled(bool(valid and self.preview.validated_revision == revision and self.preview.compiled is not None))

    def set_parameter(self, key, value):
        if self._syncing:
            return
        document = deepcopy(self.document)
        if key == "period_um":
            scale = value / document.period_um
            document.points = [[x * scale, z] for x, z in document.points]
            for region in getattr(document, "regions", []):
                region["vertices"] = [[x * scale, z] for x, z in region["vertices"]]
        setattr(document, key, value)
        self._set_document(document, fit=key in ("period_um", "d1_um", "d2_um", "zeta_deg"))

    def set_one_d_mode(self, mode):
        if self._syncing or self.document.dimension != 1 or mode == getattr(self.document, "one_d_mode", "surface"):
            return
        from custom_profile import compile_profile, surface_to_regions

        document = deepcopy(self.document)
        try:
            if mode == "regions" and not document.regions:
                document = surface_to_regions(document)
            elif mode == "native":
                compiled = compile_profile(document)
                document.native_source = compiled.source_text
                document.native_pstring = compiled.node["params"].get("pstring", "(1)")
                document.native_base_directory = ""
            elif mode == "surface" and len(document.points) < 2:
                document.points = default_profile(1).points
                ratio = document.period_um
                document.points = [[x * ratio, z] for x, z in document.points]
            document.version = 2
            document.one_d_mode = mode
        except (ValueError, OSError) as error:
            self.status_label.setText("Could not change editing mode: " + str(error))
            self._sync_controls()
            return
        self.canvas.selected_index = 0
        self.canvas.selected_vertex = 0
        self._set_document(document, fit=mode != "native")

    def load_spacer_preset(self):
        if self.document.dimension != 1:
            return
        from custom_profile import region_profile

        document = region_profile("core_spacers", medium_i=self.document.medium_i, medium_t=self.document.medium_t)
        document.nlayers = self.document.nlayers
        self.canvas.selected_index = 0
        self.canvas.selected_vertex = 0
        self._set_document(document, fit=True)

    def set_region_parameter(self, key, value):
        if self._syncing or not region_mode(self.document) or self.canvas.selected_index < 0:
            return
        document = deepcopy(self.document)
        document.regions[self.canvas.selected_index][key] = value
        self._set_document(document)

    def set_region_transform(self, key, value):
        if self._syncing or not region_mode(self.document) or self.canvas.selected_index < 0:
            return
        document = deepcopy(self.document)
        region = document.regions[self.canvas.selected_index]
        base = region_base_layer(document, region)
        if base is not None and key == "center_z":
            return
        if "base_layer" in region:
            if key != "height" or value <= 0:
                return
            region["base_layer"]["thickness"] = float(value)
            self._set_document(document)
            return
        vertices = np.asarray(region["vertices"], dtype=float)
        center = (vertices.min(axis=0) + vertices.max(axis=0)) / 2
        if key in ("center_x", "center_z"):
            axis = 0 if key == "center_x" else 1
            vertices[:, axis] += value - center[axis]
        else:
            axis = 0 if key == "width" else 1
            extent = np.ptp(vertices[:, axis])
            if extent <= 0 or value <= 0:
                return
            pivot = float(base["base_layer"]["thickness"]) if base is not None and key == "height" else center[axis]
            vertices[:, axis] = pivot + (vertices[:, axis] - pivot) * value / extent
        region["vertices"] = vertices.tolist()
        if base is not None and key == "height":
            for attachment in region.get("attachments", []):
                if "edge" in attachment and attachment["target"] != base["id"]:
                    enforce_region_constraints(document, previous=self.document,
                                               edited=(region["id"], attachment["vertex"]), axis="z")
        self._set_document(document)

    def rotate_region(self, degrees):
        if not region_mode(self.document) or self.canvas.selected_index < 0 or not degrees:
            return
        document = deepcopy(self.document)
        region = document.regions[self.canvas.selected_index]
        if "base_layer" in region:
            self.status_label.setText("The base film stays aligned with the substrate and period boundaries.")
            return
        vertices = np.asarray(region["vertices"], dtype=float)
        center = (vertices.min(axis=0) + vertices.max(axis=0)) / 2
        angle = math.radians(degrees)
        rotation = np.array([[math.cos(angle), math.sin(angle)], [-math.sin(angle), math.cos(angle)]])
        region["vertices"] = ((vertices - center) @ rotation + center).tolist()
        self._set_document(document)

    def set_region_vertex(self, axis, value):
        if self._syncing or not region_mode(self.document) or self.canvas.selected_index < 0 or self.canvas.selected_vertex < 0:
            return
        document = deepcopy(self.document)
        index, vertex = self.canvas.selected_index, self.canvas.selected_vertex
        x, z = document.regions[index]["vertices"][vertex]
        move_region_vertex(document, index, vertex, value if axis == "x" else x, value if axis == "z" else z, self.canvas.link_shared_vertices, axis=axis)
        self._set_document(document)

    def add_region(self, shape="rectangle"):
        if not region_mode(self.document):
            return
        period = self.document.period_um
        width, height = period * .2, period * .15
        center = np.array([0., period * .25])
        if shape == "film":
            vertices = [[-period / 2, 0], [period / 2, 0], [period / 2, self.document.film_thickness_um], [-period / 2, self.document.film_thickness_um]]
        elif shape == "ellipse":
            angles = np.arange(self.document.curve_segments) * (2 * math.pi / self.document.curve_segments)
            vertices = (np.column_stack((np.cos(angles) * width / 2, np.sin(angles) * height / 2)) + center).tolist()
        else:
            top_width = width * .6 if shape == "trapezoid" else width
            vertices = (np.array([[-width / 2, -height / 2], [width / 2, -height / 2], [top_width / 2, height / 2], [-top_width / 2, height / 2]]) + center).tolist()
        self.add_drawn_region(vertices, "Film" if shape == "film" else shape.title())

    def add_drawn_region(self, vertices, name="Polygon"):
        if not region_mode(self.document) or len(vertices) < 3:
            return
        document = deepcopy(self.document)
        document.regions.append({"id": str(uuid.uuid4()), "name": name + " " + str(len(document.regions) + 1), "material": document.film_material, "vertices": deepcopy(vertices)})
        self.canvas.selected_index = len(document.regions) - 1
        self.canvas.selected_vertex = 0
        self._set_document(document)

    def begin_region_drawing(self):
        self.canvas.begin_polygon()
        self.instructions.setText("Click polygon vertices in the drawing. Right-click, double-click, or press Enter to finish; Escape cancels. At least three vertices are required.")

    def duplicate_region(self):
        index = self.canvas.selected_index
        if not region_mode(self.document) or index < 0:
            return
        document = deepcopy(self.document)
        region = deepcopy(document.regions[index])
        region.pop("attachments", None)
        region.pop("base_layer", None)
        region["id"] = str(uuid.uuid4())
        region["name"] = region.get("name", "Region") + " copy"
        region["vertices"] = (np.asarray(region["vertices"]) + [document.period_um * .1, 0]).tolist()
        document.regions.append(region)
        self.canvas.selected_index = len(document.regions) - 1
        self._set_document(document)

    def delete_region(self):
        index = self.canvas.selected_index
        if not region_mode(self.document) or index < 0:
            return
        reason = region_delete_error(self.document, index)
        if reason:
            self.status_label.setText(reason)
            return
        document = deepcopy(self.document)
        del document.regions[index]
        self.canvas.selected_index = min(index, len(document.regions) - 1)
        self._set_document(document)

    def reorder_region(self, offset):
        index = self.canvas.selected_index
        if not region_mode(self.document) or not 0 <= index + offset < len(self.document.regions):
            return
        document = deepcopy(self.document)
        document.regions.insert(index + offset, document.regions.pop(index))
        self.canvas.selected_index = index + offset
        self._set_document(document)

    def add_region_vertex(self):
        index = self.canvas.selected_index
        if not region_mode(self.document) or index < 0:
            return
        document = deepcopy(self.document)
        vertices = document.regions[index]["vertices"]
        vertex = max(self.canvas.selected_vertex, 0)
        reason = region_vertex_edit_error(document, index, vertex, "insert")
        if reason:
            self.status_label.setText(reason)
            return
        a, b = np.asarray(vertices[vertex]), np.asarray(vertices[(vertex + 1) % len(vertices)])
        vertices.insert(vertex + 1, ((a + b) / 2).tolist())
        remap_region_vertices(document, document.regions[index]["id"], vertex + 1, 1)
        self.canvas.selected_vertex = vertex + 1
        self._set_document(document)

    def delete_region_vertex(self):
        index = self.canvas.selected_index
        if not region_mode(self.document) or index < 0 or len(self.document.regions[index]["vertices"]) <= 3:
            return
        document = deepcopy(self.document)
        vertex = max(self.canvas.selected_vertex, 0)
        reason = region_vertex_edit_error(document, index, vertex, "delete")
        if reason:
            self.status_label.setText(reason)
            return
        del document.regions[index]["vertices"][vertex]
        remap_region_vertices(document, document.regions[index]["id"], vertex, -1)
        self.canvas.selected_vertex = 0
        self._set_document(document)

    def set_point_coordinate(self, key, value):
        if self._syncing or self.document.dimension != 1 or self.canvas.selected_index < 0:
            return
        document = deepcopy(self.document)
        index = self.canvas.selected_index
        x, z = document.points[index]
        move_profile_point(document, index, value if key == "point_x" else x, value if key == "point_z" else z)
        self._set_document(document)

    def set_feature_parameter(self, key, value):
        if self._syncing or self.document.dimension != 2 or self.canvas.selected_index < 0:
            return
        document = deepcopy(self.document)
        document.features[self.canvas.selected_index][key] = value
        self._set_document(document)

    def set_vertex_coordinate(self, key, value):
        if self._syncing or self.canvas.selected_index < 0:
            return
        document = deepcopy(self.document)
        feature = document.features[self.canvas.selected_index]
        if feature["shape"] != "polygon" or self.canvas.selected_vertex < 0:
            return
        feature["vertices"][self.canvas.selected_vertex][0 if key == "vertex_x" else 1] = value
        self._set_document(document)

    def add_point(self):
        if self.document.dimension != 1:
            return
        document = deepcopy(self.document)
        index = min(max(self.canvas.selected_index, 0), len(document.points) - 2)
        a, b = document.points[index:index + 2]
        document.points.insert(index + 1, [(a[0] + b[0]) / 2, (a[1] + b[1]) / 2])
        self.canvas.selected_index = index + 1
        self._set_document(document)

    def delete_point(self):
        index = self.canvas.selected_index
        if self.document.dimension != 1 or not 0 < index < len(self.document.points) - 1:
            return
        document = deepcopy(self.document)
        del document.points[index]
        self.canvas.selected_index = max(index - 1, 0)
        self._set_document(document)

    def add_feature(self, shape):
        if self.document.dimension != 2:
            return
        document = deepcopy(self.document)
        size = min(document.d1_um, document.d2_um) * .12
        feature = {"id": str(uuid.uuid4()), "name": f"{shape.title()} {len(document.features) + 1}", "shape": shape, "role": "hole", "x_um": 0., "y_um": 0., "rotation_deg": 0., "radius_um": size, "width_um": size * 2, "height_um": size * 2, "vertices": [[-size, -size], [size, -size], [size, size], [-size, size]], "depth_um": document.film_thickness_um, "post_height_um": .1, "material": document.film_material}
        if document.features:
            feature["x_um"] = document.d1_um * .34
        document.features.append(feature)
        self.canvas.selected_index = len(document.features) - 1
        self.canvas.selected_vertex = 0
        self._set_document(document)

    def duplicate_feature(self):
        index = self.canvas.selected_index
        if self.document.dimension != 2 or index < 0:
            return
        document = deepcopy(self.document)
        feature = deepcopy(document.features[index])
        feature["id"] = str(uuid.uuid4())
        feature["name"] = feature.get("name", feature["shape"].title()) + " copy"
        feature["x_um"] += min(document.d1_um, document.d2_um) * .2
        document.features.append(feature)
        self.canvas.selected_index = len(document.features) - 1
        self._set_document(document)

    def delete_feature(self):
        index = self.canvas.selected_index
        if self.document.dimension != 2 or index < 0:
            return
        document = deepcopy(self.document)
        del document.features[index]
        self.canvas.selected_index = min(index, len(document.features) - 1)
        self._set_document(document)

    def add_vertex(self):
        document = deepcopy(self.document)
        if document.dimension != 2 or self.canvas.selected_index < 0:
            return
        feature = document.features[self.canvas.selected_index]
        if feature["shape"] != "polygon":
            return
        index = max(self.canvas.selected_vertex, 0)
        a, b = feature["vertices"][index], feature["vertices"][(index + 1) % len(feature["vertices"])]
        feature["vertices"].insert(index + 1, [(a[0] + b[0]) / 2, (a[1] + b[1]) / 2])
        self.canvas.selected_vertex = index + 1
        self._set_document(document)

    def delete_vertex(self):
        document = deepcopy(self.document)
        if document.dimension != 2 or self.canvas.selected_index < 0:
            return
        feature = document.features[self.canvas.selected_index]
        if feature["shape"] != "polygon" or len(feature["vertices"]) <= 3:
            return
        del feature["vertices"][max(self.canvas.selected_vertex, 0)]
        self.canvas.selected_vertex = 0
        self._set_document(document)

    def _begin_drag(self):
        self._drag_before = deepcopy(self.document)

    def _drag_changed(self, document):
        self._set_document(document, record=False)

    def _finish_drag(self):
        if self._drag_before is not None and document_signature(self._drag_before) != document_signature(self.document):
            self._history = self._history[:self._history_index + 1] + [deepcopy(self.document)]
            self._history_index += 1
            self._sync_controls()
        self._drag_before = None

    def undo(self):
        if self._history_index > 0:
            self._history_index -= 1
            self.document = deepcopy(self._history[self._history_index])
            self._sync_controls()
            self._update_preview()

    def redo(self):
        if self._history_index < len(self._history) - 1:
            self._history_index += 1
            self.document = deepcopy(self._history[self._history_index])
            self._sync_controls()
            self._update_preview()

    def reset(self):
        document = default_editor_profile(self.document.dimension, medium_i=self.document.medium_i, medium_t=self.document.medium_t)
        self.canvas.selected_index = 0
        self._set_document(document, fit=True)

    def reload_document(self, document):
        document = deepcopy(document)
        ensure_spacer_constraints(document)
        enforce_region_constraints(document)
        self.document = deepcopy(document)
        self._applied_document = deepcopy(document)
        self._applied_signature = document_signature(document)
        self._history = [deepcopy(document)]
        self._history_index = 0
        self.canvas.selected_index = 0
        self.canvas.selected_vertex = -1
        self._sync_controls()
        self.canvas.fit()
        self._update_preview()

    def apply(self):
        if not self.apply_button.isEnabled() or self.preview.compiled is None:
            return False
        if not self.preview.is_validated():
            self._update_preview()
            return False
        if document_signature(self.document) != document_signature(self.preview.document):
            self.apply_button.setEnabled(False)
            self.status_label.setText("The draft changed. Wait for the current geometry to finish validation.")
            return False
        try:
            if self.apply_callback is not None:
                self.apply_callback(deepcopy(self.document), deepcopy(self.preview.compiled.node))
        except Exception as exc:
            self.status_label.setText(f"Could not apply geometry: {exc}")
            return False
        self._applied_document = deepcopy(self.document)
        self._applied_signature = document_signature(self.document)
        self._shown = True
        self._sync_controls()
        self.status_label.setText("Applied to simulation inputs. The simulation has not been started.")
        return True

    def _adapt_layout(self):
        if not hasattr(self, "preview"):
            return
        compact = self.width() < 1280
        if compact == self._compact_layout:
            return
        self._compact_layout = compact
        if compact:
            self.visual_tabs.addTab(self._drawing_panel, "Edit profile")
            self.visual_tabs.addTab(self.preview, "3D preview")
            self.visual_stack.setCurrentWidget(self.visual_tabs)
        else:
            while self.visual_tabs.count():
                self.visual_tabs.removeTab(0)
            self.visual_splitter.addWidget(self._drawing_panel)
            self.visual_splitter.addWidget(self.preview)
            self.visual_splitter.setSizes([510, 650])
            self.visual_stack.setCurrentWidget(self.visual_splitter)
            self._drawing_panel.show()
            self.preview.show()

    def resizeEvent(self, event):
        super().resizeEvent(event)
        self._adapt_layout()

    def showEvent(self, event):
        super().showEvent(event)
        self._shown = True
        self._update_preview()

    def closeEvent(self, event):
        if self.dirty:
            choice = QMessageBox.question(self, "Unapplied profile changes", "Discard the changes that have not been applied to the simulation?", QMessageBox.Discard | QMessageBox.Cancel, QMessageBox.Cancel)
            if choice != QMessageBox.Discard:
                event.ignore()
                return
            self.reload_document(self._applied_document)
        self.preview.cancel()
        event.accept()

    def shutdown(self):
        self.canvas._draw_pending = False
        self.preview.shutdown()
