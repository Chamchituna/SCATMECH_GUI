from __future__ import annotations
import csv
import operator
import numpy as np


def _read_numeric_csv(path: str) -> np.ndarray:

    rows = []
    with open(path, "r", encoding="utf-8-sig") as f:
        rdr = csv.reader(f)
        for line_number, raw in enumerate(rdr, start=1):
            if not raw or not any(cell.strip() for cell in raw):
                continue
            parts = raw if len(raw) > 1 else raw[0].split()
            try:
                values = [float(x) for x in parts]
            except ValueError as exc:
                if not rows and [part.strip().lower() for part in parts] == ["theta", "rp", "rs"]:
                    continue
                raise ValueError(f"Nonnumeric reflectance data at line {line_number}.") from exc
            if rows and len(values) != len(rows[0]):
                raise ValueError(f"Inconsistent number of columns at line {line_number}.")
            rows.append(values)
    if not rows:
        raise ValueError("No numeric rows found in file.")
    return np.array(rows, dtype=float)


def plot_reflectance(ax, csv_path: str, component: str = "p", semilogy: bool = False):
    if component.lower().startswith("p"):
        y_col = 1
        label = "Rp (p-polarized)"
    elif component.lower().startswith("s"):
        y_col = 2
        label = "Rs (s-polarized)"
    else:
        raise ValueError("Reflectance component must be 'p' or 's'.")
    plot_csv(ax, csv_path, y_col=y_col, semilogy=semilogy, label=label)
    ax.set_ylabel("Reflectance")

def plot_csv(ax, csv_path: str, x_col: int = 0, y_col: int = 1, semilogy: bool = False, label: str | None = None):

    data = _read_numeric_csv(csv_path)
    ncols = data.shape[1]
    try:
        x_col, y_col = operator.index(x_col), operator.index(y_col)
    except TypeError as exc:
        raise ValueError("Plot columns must be nonnegative integers.") from exc
    if min(x_col, y_col) < 0:
        raise ValueError("Plot columns must be nonnegative integers.")
    if max(x_col, y_col) >= ncols:
        raise ValueError(f"CSV has {ncols} columns; need at least {max(x_col, y_col)+1}.")

    x = data[:, x_col]
    y = data[:, y_col]
    mask = np.isfinite(x) & np.isfinite(y)
    if semilogy:
        mask &= y > 0
    if not np.any(mask):
        raise ValueError("No finite plottable values in the selected columns.")
    x, y = x[mask], y[mask]

    ax.cla()
    if semilogy:
        ax.semilogy(x, y, label=label)
    else:
        ax.plot(x, y, label=label)

    ax.set_xlabel("Incidence angle θ (deg)" if x_col == 0 else "X")
    if label in ("Rp", "Rs"):
        ax.set_ylabel("Reflectance")
    else:
        ax.set_ylabel("Value")

    if label:
        ax.legend(loc="best")
    ax.grid(True, alpha=0.3)
