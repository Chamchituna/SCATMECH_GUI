from __future__ import annotations

from collections import Counter
import math
import re
from pathlib import Path
from typing import Optional, Tuple

import numpy as np
from matplotlib.ticker import AutoMinorLocator

OUTPUT_COLUMNS = ("theta_s", "phi_s", "eta", "DOLP", "DOCP", "DOP", "I")
PLOT_QUANTITIES = (
    ("auto", "Recommended for model"),
    ("intensity", "Scattering intensity"),
    ("alpha_beta", "α and β — normalized Stokes"),
    ("polarization", "DOLP, DOCP, and DOP"),
    ("dolp", "Linear polarization — DOLP"),
    ("docp", "Circular polarization — DOCP"),
    ("dop", "Total polarization — DOP"),
    ("eta", "Polarization angle η"),
)


def recommended_quantity(meta=None):
    model = (meta or {}).get("model_name") or (meta or {}).get("model_family")
    return "alpha_beta" if model in ("RCW_BRDF_Model", "CrossRCW_BRDF_Model") else "intensity"


def _scattering_name(meta):
    code = str((meta or {}).get("direction_code", ""))
    transmission = code in ("1", "3") or not code and "transmission" in str((meta or {}).get("direction_label", "")).lower()
    return "BTDF" if transmission else "BRDF"


def _incident_label(meta):
    polarization = (meta or {}).get("incident_polarization", {"label": "p"})
    label = polarization.get("label", "p") if isinstance(polarization, dict) else str(polarization)
    return label + "-polarized incidence" if label in ("p", "s") else str(label) + " incidence"


def quantity_description(quantity="auto", meta=None):
    key = recommended_quantity(meta) if quantity in (None, "auto") else quantity
    incidence = _incident_label(meta)
    if key == "alpha_beta":
        return f"α = Q/I and β = U/I for {incidence}; dimensionless, shown only where I > 0. The model’s cone half-angle α is a separate input."
    if key == "intensity":
        return f"{_scattering_name(meta)} for {incidence} [sr⁻¹], versus the sampled scattering angle."
    if key == "eta":
        return f"Polarization ellipse orientation η [degrees] for {incidence}; undefined without a linearly polarized signal."
    return f"Polarization fractions for {incidence}: DOLP and DOP range from 0 to 1; signed DOCP ranges from −1 to 1. Undefined where I ≤ 0."


def _polarization_values(data):
    valid = np.isfinite(data[:, 6]) & (data[:, 6] > 0)
    eta, dolp, docp, dop = (data[:, index].copy() for index in (2, 3, 4, 5))
    with np.errstate(invalid="ignore"):
        phase = np.deg2rad(2 * eta)
        alpha, beta = dolp * np.cos(phase), dolp * np.sin(phase)
    alpha[dolp == 0] = 0
    beta[dolp == 0] = 0
    eta[dolp <= 0] = np.nan
    values = {"alpha": alpha, "beta": beta, "eta": eta, "dolp": dolp, "docp": docp, "dop": dop}
    for value in values.values():
        value[~valid | ~np.isfinite(value)] = np.nan
    return values


def _load_numeric_table(csv_path: str) -> np.ndarray:
    rows = []
    with Path(csv_path).open("r", encoding="utf-8") as handle:
        for raw in handle:
            if not raw.strip() or raw.lstrip().startswith("#"):
                continue
            tokens = [token for token in re.split(r"[\s,]+", raw.strip()) if token]
            try:
                row = [float(token.replace("D", "E").replace("d", "e")) for token in tokens]
            except ValueError:
                for token in tokens:
                    try:
                        float(token.replace("D", "E").replace("d", "e"))
                    except ValueError:
                        continue
                    raise ValueError("BRDF output contains a malformed numeric row.") from None
                continue
            if len(row) >= 2:
                rows.append(row)
    if not rows:
        raise ValueError("No numeric rows detected in BRDF output.")
    widths = Counter(map(len, rows))
    if len(widths) != 1:
        raise ValueError("BRDF output contains rows with inconsistent column counts.")
    return np.asarray(rows, dtype=float)


def plot_csv(
    ax,
    csv_path: str,
    *,
    x_col: Optional[int] = None,
    y_col: Optional[int] = None,
    semilogy: Optional[bool] = None,
    x_span_hint: Tuple[float, float] = (-90.0, 90.0),
    label: Optional[str] = None,
    meta: Optional[dict] = None,
    quantity: str = "auto",
) -> None:
    data = _load_numeric_table(csv_path)
    ncols = data.shape[1]
    documented = ncols == len(OUTPUT_COLUMNS)
    names = list(OUTPUT_COLUMNS) if documented else list((meta or {}).get("header_tokens") or [])
    names += [f"C{i + 1}" for i in range(len(names), ncols)]
    for column, name in ((x_col, "X"), (y_col, "Y")):
        if column is not None and (not isinstance(column, (int, np.integer)) or not 0 <= column < ncols):
            raise ValueError(f"{name} column must be between 0 and {ncols - 1}.")
    if x_col is None:
        x_col = 0
        if documented:
            theta = data[:, 0][np.isfinite(data[:, 0])]
            phi = data[:, 1][np.isfinite(data[:, 1])]
            if theta.size and phi.size and np.ptp(theta) < 1e-9 and np.ptp(phi) > 1e-9:
                x_col = 1
    if x_col == y_col:
        raise ValueError("X and Y columns must be different.")
    mask = np.isfinite(data[:, x_col])
    group_col = (1 - x_col) if documented and x_col in (0, 1) else None
    if group_col is not None:
        mask &= np.isfinite(data[:, group_col])
    data = data[mask]
    if not len(data):
        raise ValueError("No finite scan coordinates available for plotting.")
    key = recommended_quantity(meta) if quantity in (None, "auto") else quantity
    if key not in dict(PLOT_QUANTITIES):
        raise ValueError("Unknown BRDF plot quantity: " + str(key))
    polarization = _polarization_values(data) if documented else {}
    if y_col is not None:
        explicit = {2: "eta", 3: "dolp", 4: "docp", 5: "dop", 6: "intensity"}
        key = explicit.get(y_col, "column") if documented else "column"
        values = polarization[key] if key in polarization else data[:, y_col]
        series = [(names[y_col], values, "tab:blue", "-")]
    elif key == "intensity":
        y_col = ncols - 1
        if y_col == x_col:
            raise ValueError("X and Y columns must be different.")
        series = [(names[y_col], data[:, y_col], "tab:blue", "-")]
    elif not documented:
        raise ValueError("Polarization plots require the seven-column BRDF output table.")
    elif key == "alpha_beta":
        series = [(r"$\alpha$ = Q/I", polarization["alpha"], "tab:orange", "-"),
                  (r"$\beta$ = U/I", polarization["beta"], "tab:blue", "--")]
    else:
        chosen = ("dolp", "docp", "dop") if key == "polarization" else (key,)
        series = [(r"$\eta$" if name == "eta" else name.upper(), polarization[name], color, style)
                  for name, color, style in zip(chosen, ("tab:blue", "tab:orange", "tab:green"), ("-", "--", ":"))]
    groups = np.unique(data[:, group_col]) if group_col is not None else np.array([0.0])
    finite_values = np.concatenate([values[np.isfinite(values)] for name, values, color, style in series])
    normalized = key in ("alpha_beta", "polarization", "dolp", "docp", "dop")
    use_log = False
    if key in ("intensity", "column") and finite_values.size and np.all(finite_values > 0):
        use_log = bool(semilogy) if semilogy is not None else float(np.max(finite_values) / np.min(finite_values)) >= 100
    ax.clear()
    ax.set_yscale("log" if use_log else "linear")
    for name, values, color, style in series:
        for group_index, group in enumerate(groups):
            positions = np.flatnonzero(data[:, group_col] == group) if group_col is not None else np.arange(len(data))
            positions = positions[np.argsort(data[positions, x_col], kind="stable")]
            x = data[positions, x_col]
            y = np.where(np.isfinite(values[positions]), values[positions], np.nan)
            finite = np.isfinite(y)
            if not np.any(finite):
                continue
            trace_label = label or name
            if label and len(series) > 1:
                trace_label += " · " + name
            if len(groups) > 1:
                symbol = "φ" if group_col == 1 else "θ"
                trace_label += f" ({symbol} = {group:g}°)"
            isolated = np.flatnonzero(finite & ~np.r_[False, finite[:-1]] & ~np.r_[finite[1:], False])
            line_color = color if len(groups) == 1 else f"C{group_index % 10}"
            ax.plot(x, y, linewidth=1.7, color=line_color, linestyle=style,
                    marker="o" if len(isolated) else None, markevery=isolated.tolist(), label=trace_label)
    if documented and x_col in (0, 1):
        ax.set_xlabel(r"$\theta_s$ (degrees)" if x_col == 0 else r"$\phi_s$ (degrees)")
    else:
        ax.set_xlabel(names[x_col])
    if documented and key == "intensity":
        ylabel = _scattering_name(meta) + r" [sr$^{-1}$]"
    elif key == "alpha_beta":
        ylabel = "Normalized Stokes parameter (dimensionless)"
    elif normalized:
        ylabel = "Polarization fraction (dimensionless)"
    elif key == "eta":
        ylabel = r"Polarization angle $\eta$ (degrees)"
    elif documented and y_col in (0, 1):
        ylabel = r"$\theta_s$ (degrees)" if y_col == 0 else r"$\phi_s$ (degrees)"
    else:
        ylabel = names[y_col]
    ax.set_ylabel(ylabel)
    if ax.lines and (len(groups) > 1 or label or documented and key != "intensity"):
        ax.legend(loc="best")
    if normalized:
        low = 0 if key in ("dolp", "dop") else -1
        high = 1
        if finite_values.size:
            low, high = min(low, float(finite_values.min())), max(high, float(finite_values.max()))
        ax.set_ylim(low - .05, high + .05)
    if not ax.lines:
        message = "No nonzero scattered signal at the sampled angles.\nNormalized polarization is undefined; select Scattering intensity to inspect I."
        if key == "eta" and np.any(data[:, 6] > 0):
            message = "No linearly polarized signal at the sampled angles.\nPolarization orientation η is undefined."
        if key in ("intensity", "column"):
            raise ValueError("No finite angle and measurement pairs available for plotting.")
        ax.text(.5, .5, message, ha="center", va="center", wrap=True, transform=ax.transAxes)
    if documented:
        model = (meta or {}).get("model_name") or ""
        title = "α and β" if key == "alpha_beta" else _scattering_name(meta) if key == "intensity" else "Polarization"
        ax.set_title((str(model) + " · " if model else "") + title + "\n" + _incident_label(meta), fontsize=10)

    lo = float(np.min(data[:, x_col]))
    hi = float(np.max(data[:, x_col]))

    scatter = (meta or {}).get("scatter") or {}
    prefix = "azimuth" if documented and x_col == 1 else "scatter"
    try:
        hint_lo = float(scatter.get(f"{prefix}_start", lo))
        hint_hi = float(scatter.get(f"{prefix}_end", hi))
        if math.isfinite(hint_lo) and math.isfinite(hint_hi) and hint_lo < hint_hi and hint_lo <= lo and hint_hi >= hi:
            lo, hi = hint_lo, hint_hi
    except (TypeError, ValueError):
        pass
    if lo == hi:
        margin = max(abs(lo) * 0.02, 0.5)
        lo, hi = lo - margin, hi + margin
    ax.set_xlim(lo, hi)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    if not use_log:
        ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.grid(True, which="major", linewidth=0.8, alpha=0.6)
    ax.grid(True, which="minor", linewidth=0.4, alpha=0.35)
    ax.set_axisbelow(True)
    ax.figure.tight_layout()
