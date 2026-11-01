from pathlib import Path
from copy import deepcopy

import numpy as np
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QFontDatabase
from PyQt5.QtWidgets import (
    QCheckBox, QComboBox, QDoubleSpinBox, QFileDialog, QFormLayout,
    QGroupBox, QHBoxLayout, QLabel, QLineEdit, QPlainTextEdit, QPushButton,
    QSizePolicy, QVBoxLayout, QWidget,
)

from custom_profile_constraints import region_vertex_edit_error


def region_base_layer(document, region):
    regions = {item["id"]: item for item in document.regions}
    pending = [region]
    visited = set()
    while pending:
        current = pending.pop()
        identity = current["id"]
        if identity in visited:
            continue
        visited.add(identity)
        if "base_layer" in current:
            return current
        pending.extend(regions[item["target"]] for item in current.get("attachments", []) if item.get("target") in regions)
    return None


class RegionProperties(QGroupBox):
    def __init__(self, editor):
        super().__init__("Selected material region")
        self.editor = editor
        form = QFormLayout(self)
        form.setRowWrapPolicy(QFormLayout.WrapLongRows)
        self.selection = QComboBox()
        self.selection.setSizePolicy(QSizePolicy.Ignored, QSizePolicy.Fixed)
        self.selection.setMinimumWidth(80)
        self.selection.currentIndexChanged.connect(lambda index: editor.canvas.set_selection(index, 0) if not editor._syncing else None)
        form.addRow("Region", self.selection)
        self.name = QLineEdit()
        self.name.setSizePolicy(QSizePolicy.Ignored, QSizePolicy.Fixed)
        self.name.setMinimumWidth(50)
        self.name.editingFinished.connect(lambda: editor.set_region_parameter("name", self.name.text()))
        form.addRow("Name", self.name)
        self.material_kind = QComboBox()
        self.material_kind.setSizePolicy(QSizePolicy.Ignored, QSizePolicy.Fixed)
        self.material_kind.setMinimumWidth(90)
        self.material_kind.addItems(["Optical material", "Incident medium / void", "Substrate material"])
        self.material_kind.currentIndexChanged.connect(self._material_kind)
        form.addRow("Fill", self.material_kind)
        self.material = editor._material(form, "Material", lambda value: editor.set_region_parameter("material", value))
        self.coordinates = {}
        for key, label in (("center_x", "Center x (µm)"), ("center_z", "Center z (µm)"), ("width", "Width (µm)"), ("height", "Height (µm)")):
            field = self._number()
            if key in ("width", "height"):
                field.setMinimum(.00000001)
            field.valueChanged.connect(lambda value, name=key: editor.set_region_transform(name, value))
            self.coordinates[key] = field
            form.addRow(label, field)
        rotation = QWidget()
        rotation_row = QHBoxLayout(rotation)
        rotation_row.setContentsMargins(0, 0, 0, 0)
        self.rotation = self._number()
        self.rotation.setRange(-360, 360)
        self.rotate_button = QPushButton("Rotate")
        self.rotate_button.clicked.connect(lambda: editor.rotate_region(self.rotation.value()))
        rotation_row.addWidget(self.rotation)
        rotation_row.addWidget(self.rotate_button)
        form.addRow("Angle (degrees)", rotation)
        self.vertex = QComboBox()
        self.vertex.setSizePolicy(QSizePolicy.Ignored, QSizePolicy.Fixed)
        self.vertex.setMinimumWidth(80)
        self.vertex.currentIndexChanged.connect(lambda index: editor.canvas.set_selection(editor.canvas.selected_index, index) if not editor._syncing else None)
        form.addRow("Vertex", self.vertex)
        self.vertex_x = self._number()
        self.vertex_z = self._number()
        self.vertex_x.valueChanged.connect(lambda value: editor.set_region_vertex("x", value))
        self.vertex_z.valueChanged.connect(lambda value: editor.set_region_vertex("z", value))
        form.addRow("Vertex x (µm)", self.vertex_x)
        form.addRow("Vertex z (µm)", self.vertex_z)
        self.constraint = QLabel()
        self.constraint.setTextFormat(Qt.PlainText)
        self.constraint.setWordWrap(True)
        form.addRow("Connection", self.constraint)
        actions = QWidget()
        actions_row = QHBoxLayout(actions)
        actions_row.setContentsMargins(0, 0, 0, 0)
        self.add_vertex = QPushButton("Add vertex")
        self.delete_vertex = QPushButton("Delete vertex")
        self.add_vertex.clicked.connect(editor.add_region_vertex)
        self.delete_vertex.clicked.connect(editor.delete_region_vertex)
        actions_row.addWidget(self.add_vertex)
        actions_row.addWidget(self.delete_vertex)
        form.addRow(actions)
        ordering = QWidget()
        ordering_row = QHBoxLayout(ordering)
        ordering_row.setContentsMargins(0, 0, 0, 0)
        self.backward = QPushButton("Send backward")
        self.forward = QPushButton("Bring forward")
        self.backward.clicked.connect(lambda: editor.reorder_region(-1))
        self.forward.clicked.connect(lambda: editor.reorder_region(1))
        ordering_row.addWidget(self.backward)
        ordering_row.addWidget(self.forward)
        form.addRow(ordering)
        self.snap = QCheckBox("Snap to vertices and cell edges")
        self.snap.setChecked(True)
        self.snap.toggled.connect(lambda value: setattr(editor.canvas, "snap_vertices", value))
        self.link = QCheckBox("Move shared vertices together")
        self.link.setToolTip("Links other coincident vertices. Required preset connections to the core, film, and substrate remain attached when this is off.")
        self.link.setChecked(True)
        self.link.toggled.connect(lambda value: setattr(editor.canvas, "link_shared_vertices", value))
        form.addRow(self.snap)
        form.addRow(self.link)
        note = QLabel("Later regions cover earlier ones. Regions repeat across the period; negative z can cut into the substrate. Use a void region for holes and undercuts. Required contacts limit transforms; dimension fields show the resulting shape.")
        note.setWordWrap(True)
        form.addRow(note)

    def _number(self):
        field = QDoubleSpinBox()
        field.setDecimals(8)
        field.setRange(-1000000, 1000000)
        field.setSingleStep(.01)
        field.setKeyboardTracking(False)
        field.setMinimumWidth(70)
        field.setSizePolicy(QSizePolicy.Ignored, QSizePolicy.Fixed)
        return field

    def _material_kind(self, index):
        if self.editor._syncing:
            return
        value = (self.editor.document.film_material, "@incident", "@substrate")[index]
        self.editor.set_region_parameter("material", value)

    def sync(self):
        document = self.editor.document
        self.selection.clear()
        self.selection.addItems([f"{index + 1}. {region.get('name', 'Region')}" for index, region in enumerate(document.regions)])
        index = min(max(self.editor.canvas.selected_index, 0), len(document.regions) - 1)
        self.editor.canvas.selected_index = index
        self.selection.setCurrentIndex(index)
        self.setEnabled(index >= 0)
        if index < 0:
            return
        region = document.regions[index]
        anchored = "base_layer" in region
        base = region_base_layer(document, region)
        vertices = np.asarray(region["vertices"], dtype=float)
        self.name.setText(region.get("name", ""))
        material = region.get("material", document.film_material)
        self.material_kind.setCurrentIndex(1 if material == "@incident" else 2 if material == "@substrate" else 0)
        self.material.setText(material)
        self.material.parentWidget().setEnabled(not material.startswith("@"))
        center = (vertices.min(axis=0) + vertices.max(axis=0)) / 2
        width, height = np.ptp(vertices, axis=0)
        for key, value in zip(("center_x", "center_z", "width", "height"), [*center, width, height]):
            self.coordinates[key].setValue(value)
            self.coordinates[key].setEnabled((not anchored or key == "height") and not (base is not None and key == "center_z"))
            self.coordinates[key].setToolTip("The base film spans one period and rests on the substrate. Edit its height to change thickness." if anchored else "This region stays connected to the film surface. Edit its height or free vertices." if base is not None else "")
        self.rotation.setEnabled(not anchored)
        self.rotate_button.setEnabled(not anchored)
        self.rotate_button.setToolTip("The base film stays aligned with the substrate." if anchored else "")
        self.vertex.clear()
        self.vertex.addItems([f"Vertex {position + 1}" for position in range(len(vertices))])
        selected_vertex = min(max(self.editor.canvas.selected_vertex, 0), len(vertices) - 1)
        self.editor.canvas.selected_vertex = selected_vertex
        self.vertex.setCurrentIndex(selected_vertex)
        self.vertex_x.setValue(vertices[selected_vertex, 0])
        self.vertex_z.setValue(vertices[selected_vertex, 1])
        attachment = next((item for item in region.get("attachments", []) if item.get("vertex") == selected_vertex), None)
        self.vertex_x.setReadOnly(anchored)
        self.vertex_z.setReadOnly(anchored and selected_vertex < 2)
        if anchored:
            self.constraint.setText("Spans one period; attached to substrate. Edit height or a top vertex z to change film thickness.")
        elif attachment is None:
            self.constraint.setText("Free vertex")
        else:
            target = next((item for item in document.regions if item.get("id") == attachment.get("target")), {})
            name = target.get("name", "target region")
            surface = "film surface" if "base_layer" in target else "sidewall"
            self.constraint.setText(f"Attached to {name} edge — slides along {surface}" if "edge" in attachment else f"Shared with {name} vertex")
            resolved = attachment
            visited = set()
            while "target_vertex" in resolved and target and target.get("id") not in visited:
                visited.add(target.get("id"))
                following = next((item for item in target.get("attachments", []) if item.get("vertex") == resolved["target_vertex"]), None)
                if following is None:
                    break
                resolved = following
                target = next((item for item in document.regions if item.get("id") == resolved.get("target")), {})
            if "edge" in resolved and target:
                a, b = np.asarray(target["vertices"], dtype=float)[resolved["edge"]]
                tolerance = max(document.period_um, 1.) * 1e-12
                self.vertex_x.setReadOnly(abs(b[0] - a[0]) <= tolerance)
                self.vertex_z.setReadOnly(abs(b[1] - a[1]) <= tolerance)
        delete_reason = region_vertex_edit_error(document, index, selected_vertex, "delete")
        insert_reason = region_vertex_edit_error(document, index, selected_vertex, "insert")
        self.delete_vertex.setEnabled(len(vertices) > 3 and not delete_reason)
        self.delete_vertex.setToolTip(delete_reason)
        self.add_vertex.setEnabled(not insert_reason)
        self.add_vertex.setToolTip(insert_reason)
        self.backward.setEnabled(index > 0)
        self.forward.setEnabled(index < len(document.regions) - 1)


class NativeDefinitionPanel(QWidget):
    def __init__(self, editor):
        super().__init__()
        self.editor = editor
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        actions = QHBoxLayout()
        self.import_button = QPushButton("Import definition…")
        self.export_button = QPushButton("Export definition…")
        self.import_button.clicked.connect(self.import_file)
        self.export_button.clicked.connect(self.export_file)
        actions.addWidget(self.import_button)
        actions.addWidget(self.export_button)
        actions.addStretch()
        layout.addLayout(actions)
        help_text = QLabel("SCATMECH sections: PARAMETERS, optional WORKING, MATERIALS, VERTICES, BOUNDARIES; terminate each with END. MATERIALS also supports ANISO, MAGNETIC, and ANISOMAGNETIC declarations. A boundary is: vertex1 vertex2 material_on_left material_on_right. Use medium_i / medium_t for inherited media. Coordinate expressions can reference parameters and working values. The vertex extent defines the period.")
        help_text.setWordWrap(True)
        layout.addWidget(help_text)
        self.source = QPlainTextEdit()
        self.source.setFont(QFontDatabase.systemFont(QFontDatabase.FixedFont))
        self.source.setLineWrapMode(QPlainTextEdit.NoWrap)
        self.source.textChanged.connect(self._source_changed)
        layout.addWidget(self.source, 1)
        parameters = QFormLayout()
        self.parameters = QLineEdit()
        self.parameters.setPlaceholderText("(1)")
        self.parameters.textChanged.connect(lambda value: editor.set_parameter("native_pstring", value))
        parameters.addRow("Parameter values", self.parameters)
        layout.addLayout(parameters)
        self.note = QLabel("Native definitions are evaluated by SCATMECH. The preview shows evaluated solver layers; Apply stays disabled until the current definition validates.")
        self.note.setWordWrap(True)
        layout.addWidget(self.note)

    def _source_changed(self):
        self.editor.set_parameter("native_source", self.source.toPlainText())

    def sync(self):
        document = self.editor.document
        if self.source.toPlainText() != document.native_source:
            self.source.setPlainText(document.native_source)
        self.parameters.setText(document.native_pstring)

    def import_path(self, path):
        source = Path(path).read_text(encoding="utf-8")
        document = deepcopy(self.editor.document)
        document.native_source = source
        document.native_base_directory = str(Path(path).resolve().parent)
        self.editor._set_document(document)

    def export_path(self, path):
        Path(path).write_text(self.editor.document.native_source, encoding="utf-8")

    def import_file(self):
        filename, _ = QFileDialog.getOpenFileName(self, "Import SCATMECH grating definition", "", "Grating definitions (*.grt *.dat *.txt);;All files (*)")
        if filename:
            try:
                self.import_path(filename)
            except (OSError, UnicodeError) as error:
                self.editor.status_label.setText("Could not import definition: " + str(error))

    def export_file(self):
        filename, _ = QFileDialog.getSaveFileName(self, "Export SCATMECH grating definition", "custom-grating.grt", "Grating definitions (*.grt);;All files (*)")
        if filename:
            try:
                self.export_path(filename)
            except OSError as error:
                self.editor.status_label.setText("Could not export definition: " + str(error))
