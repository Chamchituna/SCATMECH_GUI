from copy import deepcopy
import math

import numpy as np
from PyQt5.QtCore import pyqtSignal
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure
from matplotlib.path import Path as MplPath
from matplotlib.patches import Polygon


def feature_vertices(feature, segments=192):
    shape = feature["shape"]
    if shape == "circle":
        angles = np.linspace(0, 2 * math.pi, max(12, segments), endpoint=False)
        vertices = np.column_stack((np.cos(angles), np.sin(angles))) * feature["radius_um"]
    elif shape == "rectangle":
        vertices = np.array([[-1, -1], [1, -1], [1, 1], [-1, 1]], dtype=float)
        vertices *= np.array([feature["width_um"], feature["height_um"]]) / 2
    else:
        vertices = np.asarray(feature["vertices"], dtype=float)
    angle = math.radians(feature.get("rotation_deg", 0))
    rotation = np.array([[math.cos(angle), -math.sin(angle)], [math.sin(angle), math.cos(angle)]])
    return vertices @ rotation.T + [feature["x_um"], feature["y_um"]]


def profile_lattice(document):
    angle = math.radians(document.zeta_deg)
    return np.array([[document.d1_um, 0], [document.d2_um * math.sin(angle), document.d2_um * math.cos(angle)]])


def move_profile_point(document, index, x, z):
    points = document.points
    if index in (0, len(points) - 1):
        points[0][1] = points[-1][1] = float(z)
    else:
        points[index] = [min(max(float(x), points[index - 1][0]), points[index + 1][0]), float(z)]


class CustomProfileCanvas(FigureCanvasQTAgg):
    selection_changed = pyqtSignal(int, int)
    edit_started = pyqtSignal()
    geometry_edited = pyqtSignal(object)
    edit_finished = pyqtSignal()

    def __init__(self, parent=None):
        self.figure = Figure(figsize=(5, 5), layout="constrained")
        super().__init__(self.figure)
        self.setParent(parent)
        self.axes = self.figure.add_subplot(111)
        self.document = None
        self.selected_index = 0
        self.selected_vertex = -1
        self._drag = None
        self._limits = None
        self._handles = []
        self.mpl_connect("button_press_event", self._press)
        self.mpl_connect("motion_notify_event", self._motion)
        self.mpl_connect("button_release_event", self._release)
        self.mpl_connect("scroll_event", self._scroll)

    def set_document(self, document, fit=False):
        self.document = document
        count = len(document.points) if document.dimension == 1 else len(document.features)
        self.selected_index = min(max(self.selected_index, 0), count - 1)
        self.redraw(fit=fit)

    def set_selection(self, index, vertex=-1):
        self.selected_index = index
        self.selected_vertex = vertex
        self.selection_changed.emit(index, vertex)
        self.redraw()

    def fit(self):
        self.redraw(fit=True)

    def redraw(self, fit=False):
        if self.document is None:
            return
        limits = None if fit or self._limits is None else (self.axes.get_xlim(), self.axes.get_ylim())
        self.axes.clear()
        self._handles = []
        document = self.document
        if document.dimension == 1:
            points = np.asarray(document.points, dtype=float)
            self.axes.fill_between(points[:, 0], 0, points[:, 1], color="#e69f00", alpha=.4)
            self.axes.plot(points[:, 0], points[:, 1], "o-", color="#0072b2", markersize=6)
            if self.selected_index >= 0:
                point = points[self.selected_index]
                self.axes.scatter(*point, s=85, color="#d55e00", zorder=5)
            self._handles = [("point", i, point) for i, point in enumerate(points)]
            self.axes.axhline(0, color="#555555", linewidth=.8)
            self.axes.set_xlabel("x (µm)")
            self.axes.set_ylabel("z above substrate (µm)")
            self.axes.set_title("1D surface profile")
            width = max(abs(document.period_um), .001)
            low, high = min(0, np.min(points[:, 1])), max(.01, np.max(points[:, 1]))
            self.axes.set_xlim(-width * .6, width * .6)
            self.axes.set_ylim(low - max(width * .05, (high - low) * .15), high + max(width * .1, (high - low) * .2))
        else:
            lattice = profile_lattice(document)
            corners = np.array([[-.5, -.5], [.5, -.5], [.5, .5], [-.5, .5]]) @ lattice
            self.axes.add_patch(Polygon(corners, facecolor="#f6e6bf", edgecolor="#555555", linewidth=1.3))
            for index, feature in enumerate(document.features):
                vertices = feature_vertices(feature, document.curve_segments)
                selected = index == self.selected_index
                self.axes.add_patch(Polygon(vertices, facecolor="#b9dbe8" if feature["role"] == "hole" else "#009e73", edgecolor="#d55e00" if selected else "#0072b2", linewidth=2 if selected else 1.2, alpha=.8))
                self.axes.text(feature["x_um"], feature["y_um"], str(index + 1), ha="center", va="center", fontsize=9)
                if selected:
                    self._feature_handles(feature, vertices)
            extent = max(np.ptp(corners[:, 0]), np.ptp(corners[:, 1]), .01)
            self.axes.set_xlim(np.min(corners[:, 0]) - extent * .18, np.max(corners[:, 0]) + extent * .18)
            self.axes.set_ylim(np.min(corners[:, 1]) - extent * .18, np.max(corners[:, 1]) + extent * .18)
            self.axes.set_xlabel("x (µm)")
            self.axes.set_ylabel("y (µm)")
            self.axes.set_title("2D unit cell · holes in blue, posts in green")
        for kind, index, point in self._handles:
            if kind != "point":
                marker = "D" if kind == "rotate" else "s"
                fill = "#d55e00" if kind == "vertex" and index == self.selected_vertex else "white"
                self.axes.plot(*point, marker=marker, markersize=6, markerfacecolor=fill, markeredgecolor="#d55e00", zorder=6)
                if kind == "vertex":
                    self.axes.annotate(str(index + 1), point, xytext=(5, 5), textcoords="offset points", fontsize=8, color="#8f3900")
        self.axes.set_aspect("equal", adjustable="box")
        self.axes.grid(alpha=.2)
        if limits is not None:
            self.axes.set_xlim(*limits[0])
            self.axes.set_ylim(*limits[1])
        self._limits = (self.axes.get_xlim(), self.axes.get_ylim())
        self.draw_idle()

    def _feature_handles(self, feature, vertices):
        center = np.array([feature["x_um"], feature["y_um"]])
        if feature["shape"] == "circle":
            self._handles.append(("radius", -1, center + [feature["radius_um"], 0]))
        elif feature["shape"] == "rectangle":
            self._handles.append(("size", -1, vertices[2]))
        else:
            self._handles.extend(("vertex", index, point) for index, point in enumerate(vertices))
        if feature["shape"] != "circle":
            radius = max(np.max(np.linalg.norm(vertices - center, axis=1)), .01) * 1.35
            angle = math.radians(feature.get("rotation_deg", 0) + 90)
            handle = center + radius * np.array([math.cos(angle), math.sin(angle)])
            self.axes.plot([center[0], handle[0]], [center[1], handle[1]], color="#d55e00", linewidth=.8)
            self._handles.append(("rotate", -1, handle))

    def _press(self, event):
        if event.button != 1 or event.inaxes is not self.axes or self.document is None:
            return
        selected = None
        for kind, index, point in reversed(self._handles):
            pixel = self.axes.transData.transform(point)
            if np.linalg.norm(pixel - [event.x, event.y]) < 12:
                selected = (kind, index)
                break
        if self.document.dimension == 1:
            if selected is None:
                return
            self.selected_index = selected[1]
            self.selected_vertex = -1
        elif selected is None:
            for index in reversed(range(len(self.document.features))):
                feature = self.document.features[index]
                if MplPath(feature_vertices(feature)).contains_point([event.xdata, event.ydata]):
                    self.selected_index = index
                    self.selected_vertex = -1
                    selected = ("center", -1)
                    break
            if selected is None:
                return
        elif selected[0] == "vertex":
            self.selected_vertex = selected[1]
        self.selection_changed.emit(self.selected_index, self.selected_vertex)
        self._drag = {"kind": selected[0], "index": selected[1], "start": [event.xdata, event.ydata], "document": deepcopy(self.document)}
        self.edit_started.emit()
        self.redraw()

    def drag_to(self, x, y):
        if self._drag is None:
            return
        document = deepcopy(self._drag["document"])
        kind, index = self._drag["kind"], self._drag["index"]
        if document.dimension == 1:
            move_profile_point(document, index, x, y)
        else:
            feature = document.features[self.selected_index]
            center = np.array([feature["x_um"], feature["y_um"]])
            delta = np.array([x, y]) - center
            angle = math.radians(feature.get("rotation_deg", 0))
            local = np.array([[math.cos(angle), math.sin(angle)], [-math.sin(angle), math.cos(angle)]]) @ delta
            if kind == "center":
                change = np.array([x, y]) - self._drag["start"]
                feature["x_um"], feature["y_um"] = (center + change).tolist()
            elif kind == "radius":
                feature["radius_um"] = float(np.linalg.norm(delta))
            elif kind == "size":
                feature["width_um"], feature["height_um"] = (2 * abs(local)).tolist()
            elif kind == "rotate":
                feature["rotation_deg"] = math.degrees(math.atan2(delta[1], delta[0])) - 90
            elif kind == "vertex":
                feature["vertices"][index] = local.tolist()
        self.geometry_edited.emit(document)

    def _motion(self, event):
        if self._drag is not None and event.inaxes is self.axes and event.xdata is not None:
            self.drag_to(event.xdata, event.ydata)

    def _release(self, event):
        if self._drag is not None:
            self._drag = None
            self.edit_finished.emit()

    def _scroll(self, event):
        if event.inaxes is not self.axes:
            return
        factor = .85 if event.button == "up" else 1 / .85
        for getter, setter, value in ((self.axes.get_xlim, self.axes.set_xlim, event.xdata), (self.axes.get_ylim, self.axes.set_ylim, event.ydata)):
            lo, hi = getter()
            setter(value + (lo - value) * factor, value + (hi - value) * factor)
        self._limits = (self.axes.get_xlim(), self.axes.get_ylim())
        self.draw_idle()
