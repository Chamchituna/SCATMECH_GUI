from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple

import numpy as np


@dataclass
class GeometryRequest:
    grating: Dict[str, Any]
    wavelength_um: float
    incident_deg: float = 0.0
    rotation_deg: float = 0.0
    configuration: int = 0
    source_label: str = "RCW"
    resolution: int = 48


@dataclass
class GeometryLayer:
    z_top: float
    z_bottom: float
    u_edges: np.ndarray
    v_edges: np.ndarray
    tensors: np.ndarray
    basis: np.ndarray
    offset: np.ndarray = field(default_factory=lambda: np.zeros(2))
    label: str = ""
    continuous: bool = False
    radial: Optional[Dict[str, np.ndarray]] = None


@dataclass
class GeometryScene:
    model: str
    wavelength_um: float
    lattice: np.ndarray
    layers: List[GeometryLayer]
    medium_i: np.ndarray
    medium_t: np.ndarray
    dimension: int
    incident_deg: float = 0.0
    rotation_deg: float = 0.0
    configuration: int = 0
    files: List[str] = field(default_factory=list)
    notes: List[str] = field(default_factory=list)

    @property
    def thickness(self):
        return max((-layer.z_bottom for layer in self.layers), default=0.0)


@dataclass
class GeometryMesh:
    faces: np.ndarray
    material_ids: np.ndarray
    tensors: np.ndarray
    materials: List[Tuple[str, np.ndarray]]
    bounds: np.ndarray
    lattice_lines: List[np.ndarray]
    notes: List[str] = field(default_factory=list)
