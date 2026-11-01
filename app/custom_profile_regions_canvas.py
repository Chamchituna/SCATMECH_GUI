from copy import deepcopy

import numpy as np
from PyQt5.QtCore import pyqtSignal
from matplotlib import colormaps
from matplotlib.path import Path
from matplotlib.patches import Polygon

from custom_profile_canvas import CustomProfileCanvas
from custom_profile_constraints import enforce_region_constraints


def region_mode(document):
    return document.dimension == 1 and getattr(document, "one_d_mode", "surface") == "regions"


def move_region_vertex(document, region_index, vertex_index, x, z, linked=False, axis=None):
    previous = deepcopy(document)
    original = np.asarray(document.regions[region_index]["vertices"][vertex_index], dtype=float)
    target = [float(x), float(z)]
    if linked:
        tolerance = max(abs(document.period_um), 1.) * 1e-10
        for region in document.regions:
            for index, vertex in enumerate(region["vertices"]):
                if np.linalg.norm(np.asarray(vertex) - original) <= tolerance:
                    region["vertices"][index] = target.copy()
    else:
        document.regions[region_index]["vertices"][vertex_index] = target
    enforce_region_constraints(document, previous=previous,
                               edited=(document.regions[region_index]["id"], vertex_index), axis=axis)


class CustomProfileRegionsCanvas(CustomProfileCanvas):
    polygon_created = pyqtSignal(object)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.snap_vertices = True
        self.link_shared_vertices = True
        self._drawing = None
        self.mpl_connect("key_press_event", self._key)

    def set_document(self, document, fit=False):
        if not region_mode(document):
            self._drawing = None
            super().set_document(document, fit)
            return
        self.document = document
        self.selected_index = min(max(self.selected_index, 0), len(document.regions) - 1)
        self.redraw(fit=fit)

    def redraw(self, fit=False):
        if self.document is None or not region_mode(self.document):
            super().redraw(fit)
            return
        limits = None if fit or self._limits is None else (self.axes.get_xlim(), self.axes.get_ylim())
        self.axes.clear()
        self._handles = []
        period = self.document.period_um
        all_vertices = [np.asarray(region["vertices"], dtype=float) for region in self.document.regions]
        valid = [vertices for vertices in all_vertices if len(vertices)]
        low = min([0.] + [float(vertices[:, 1].min()) for vertices in valid])
        high = max([period * .1] + [float(vertices[:, 1].max()) for vertices in valid])
        margin = max(period * .08, (high - low) * .15)
        self.axes.axhspan(low - margin, 0, color="#d5d5d5", alpha=.65)
        self.axes.axvline(-period / 2, color="#777777", linestyle="--", linewidth=1)
        self.axes.axvline(period / 2, color="#777777", linestyle="--", linewidth=1)
        self.axes.axhline(0, color="#777777", linewidth=.7)
        for index, (region, vertices) in enumerate(zip(self.document.regions, all_vertices)):
            if len(vertices) < 3:
                continue
            material = region.get("material", "")
            color = "#dceef8" if material == "@incident" else colormaps["tab10"](((index + 1) % 10) / 9)
            if material == "@substrate":
                color = "#a4a4a4"
            selected = index == self.selected_index
            for shift in (-period, period):
                copied = vertices + [shift, 0]
                if copied[:, 0].max() >= -period * .6 and copied[:, 0].min() <= period * .6:
                    self.axes.add_patch(Polygon(copied, facecolor=color, edgecolor="#999999", alpha=.12, linewidth=.6))
            self.axes.add_patch(Polygon(vertices, facecolor=color, edgecolor="#bd3d00" if selected else "#444444", linewidth=2 if selected else .8, alpha=.85, hatch="///" if material == "@incident" else None))
            center = (vertices.min(axis=0) + vertices.max(axis=0)) / 2
            self.axes.text(*center, str(index + 1), ha="center", va="center", fontsize=9,
                           bbox={"facecolor": "white", "edgecolor": "none", "alpha": .7, "pad": 1})
            if selected:
                for vertex_index, vertex in enumerate(vertices):
                    self._handles.append(("region_vertex", vertex_index, vertex))
                    self.axes.plot(*vertex, "s", markersize=6, markerfacecolor="#bd3d00" if vertex_index == self.selected_vertex else "white", markeredgecolor="#bd3d00", zorder=10)
                    if len(vertices) <= 24 or vertex_index == self.selected_vertex:
                        self.axes.annotate(str(vertex_index + 1), vertex, xytext=(5, 5), textcoords="offset points", fontsize=8, color="#8f3900", zorder=11)
        if self._drawing:
            points = np.asarray(self._drawing)
            self.axes.plot(points[:, 0], points[:, 1], "o--", color="#a21caf", linewidth=1.5)
        self.axes.set_xlabel("x (µm)")
        self.axes.set_ylabel("z above substrate (µm)")
        self.axes.set_title("1D material regions · higher list numbers cover lower ones", fontsize=10)
        self.axes.set_xlim(min([-period * .6] + [float(vertices[:, 0].min()) - margin for vertices in valid]),
                           max([period * .6] + [float(vertices[:, 0].max()) + margin for vertices in valid]))
        self.axes.set_ylim(low - margin, high + margin)
        self.axes.set_aspect("equal", adjustable="box")
        self.axes.grid(alpha=.2)
        if limits is not None:
            self.axes.set_xlim(*limits[0])
            self.axes.set_ylim(*limits[1])
        self._limits = (self.axes.get_xlim(), self.axes.get_ylim())
        self.draw_idle()

    def begin_polygon(self):
        self._drawing = []
        self.setFocus()
        self.redraw()

    def finish_polygon(self):
        if self._drawing is not None and len(self._drawing) >= 3:
            vertices = deepcopy(self._drawing)
            self._drawing = None
            self.polygon_created.emit(vertices)
        self.redraw()

    def _key(self, event):
        if self._drawing is None:
            return
        if event.key == "escape":
            self._drawing = None
            self.redraw()
        elif event.key in ("enter", "return"):
            self.finish_polygon()

    def _press(self, event):
        if self.document is None or not region_mode(self.document):
            super()._press(event)
            return
        if event.inaxes is not self.axes:
            return
        if self._drawing is not None:
            if event.button == 3 or getattr(event, "dblclick", False):
                self.finish_polygon()
            elif event.button == 1:
                point = self._snap(event.xdata, event.ydata)
                if not self._drawing or np.linalg.norm(np.asarray(point) - self._drawing[-1]) > 1e-12:
                    self._drawing.append(list(point))
                self.redraw()
            return
        if event.button != 1:
            return
        selected = None
        for kind, index, point in reversed(self._handles):
            if np.linalg.norm(self.axes.transData.transform(point) - [event.x, event.y]) < 12:
                selected = (kind, index)
                break
        if selected is None:
            for index in reversed(range(len(self.document.regions))):
                if Path(np.asarray(self.document.regions[index]["vertices"])).contains_point([event.xdata, event.ydata]):
                    self.selected_index = index
                    self.selected_vertex = -1
                    selected = ("region_move", -1)
                    break
        elif selected[0] == "region_vertex":
            self.selected_vertex = selected[1]
        if selected is None:
            return
        self.selection_changed.emit(self.selected_index, self.selected_vertex)
        self._drag = {"kind": selected[0], "index": selected[1], "start": [event.xdata, event.ydata], "document": deepcopy(self.document)}
        self.edit_started.emit()
        self.redraw()

    def _snap(self, x, z, exclude=None):
        point = np.array([x, z], dtype=float)
        if not self.snap_vertices:
            return point
        choices = []
        reference = self._drag["document"] if self._drag is not None else self.document
        origin = np.asarray(reference.regions[exclude[0]]["vertices"][exclude[1]]) if exclude is not None else None
        for region_index, region in enumerate(self.document.regions):
            for vertex_index, vertex in enumerate(region["vertices"]):
                if exclude == (region_index, vertex_index):
                    continue
                if origin is not None and self.link_shared_vertices and np.linalg.norm(np.asarray(reference.regions[region_index]["vertices"][vertex_index]) - origin) <= max(self.document.period_um, 1.) * 1e-10:
                    continue
                choices.append(np.asarray(vertex, dtype=float))
        if choices:
            pixels = self.axes.transData.transform(np.asarray(choices))
            distances = np.linalg.norm(pixels - self.axes.transData.transform(point), axis=1)
            closest = int(np.argmin(distances))
            if distances[closest] < 9:
                return choices[closest]
        scale = self.axes.transData.transform([[0, 0], [1, 1]])
        tolerance = 6 / np.maximum(np.abs(scale[1] - scale[0]), 1)
        for boundary in (-self.document.period_um / 2, self.document.period_um / 2):
            if abs(point[0] - boundary) < tolerance[0]:
                point[0] = boundary
        if abs(point[1]) < tolerance[1]:
            point[1] = 0
        return point

    def drag_to(self, x, y):
        if self.document is None or not region_mode(self.document):
            super().drag_to(x, y)
            return
        if self._drag is None:
            return
        document = deepcopy(self._drag["document"])
        if self._drag["kind"] == "region_move":
            vertices = np.asarray(document.regions[self.selected_index]["vertices"])
            document.regions[self.selected_index]["vertices"] = (vertices + np.asarray([x, y]) - self._drag["start"]).tolist()
        else:
            x, y = self._snap(x, y, (self.selected_index, self._drag["index"]))
            move_region_vertex(document, self.selected_index, self._drag["index"], x, y, self.link_shared_vertices)
        self.geometry_edited.emit(document)
