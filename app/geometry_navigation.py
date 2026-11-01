import numpy as np
from mpl_toolkits.mplot3d import Axes3D


class ConstrainedAxes3D(Axes3D):
    geometry_bounds = None

    def constrain_pan(self):
        if self.geometry_bounds is None:
            return
        try:
            value = self.geometry_bounds() if callable(self.geometry_bounds) else self.geometry_bounds
            bounds = np.asarray(value, dtype=float)
        except Exception:
            return
        if bounds.shape != (3, 2) or not np.all(np.isfinite(bounds)):
            return
        for bounds_pair, getter, setter in zip(
            bounds,
            (self.get_xlim3d, self.get_ylim3d, self.get_zlim3d),
            (self.set_xlim3d, self.set_ylim3d, self.set_zlim3d),
        ):
            limits = np.asarray(getter(), dtype=float)
            if not np.all(np.isfinite(limits)):
                continue
            inverted = limits[0] > limits[1]
            low, high = np.sort(limits)
            span = high - low
            geometry_low, geometry_high = np.sort(bounds_pair)
            if span >= geometry_high - geometry_low:
                minimum, maximum = geometry_high - span, geometry_low
            else:
                minimum, maximum = geometry_low, geometry_high - span
            shifted_low = float(np.clip(low, minimum, maximum))
            shifted = (shifted_low, shifted_low + span)
            if inverted:
                shifted = shifted[::-1]
            if shifted != tuple(limits):
                setter(shifted, auto=None)

    def drag_pan(self, button, key, x, y):
        super().drag_pan(button, key, x, y)
        self.constrain_pan()
