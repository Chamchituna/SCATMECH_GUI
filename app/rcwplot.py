import csv
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np
from matplotlib.ticker import AutoMinorLocator, MultipleLocator


_GREEK_LABELS: Dict[str, str] = {
    "alpha": r"$\alpha$",
    "beta": r"$\beta$",
}

_GREEK_COLORS: Dict[str, str] = {
    "alpha": "tab:orange",
    "beta": "tab:blue",
}


def _nice_step(value: float) -> float:

    if not np.isfinite(value) or value <= 0:
        return 1.0

    exponent = float(np.floor(np.log10(value)))
    fraction = float(value / (10 ** exponent))
    for candidate in (1.0, 2.0, 2.5, 5.0, 10.0):
        if fraction <= candidate:
            return float(candidate * (10 ** exponent))
    return float(10 ** (exponent + 1))


def _read_csv(csv_path: str) -> Tuple[List[str], List[List[str]]]:

    path = Path(csv_path)
    if not path.exists():
        raise FileNotFoundError(f"CSV not found: {csv_path}")

    with open(path, newline="", encoding="utf-8") as handle:
        rows = [row for row in csv.reader(handle) if any(cell.strip() for cell in row)]

    if not rows:
        raise ValueError(f"CSV appears empty: {csv_path}")

    header = rows[0]
    data = rows[1:] if len(rows) > 1 else []
    if not data:
        raise ValueError("CSV contains header but no data rows.")
    if any(len(row) != len(header) for row in data):
        raise ValueError("RCW CSV contains an incomplete table row.")
    return header, data


def _column_by_hint(header: Sequence[str], hints: Iterable[str]) -> int:

    lowered = [name.lower() for name in header]
    for hint in hints:
        hint_lower = hint.lower()
        for idx, candidate in enumerate(lowered):
            if hint_lower in candidate:
                return idx
    return -1


def _to_float(values: Sequence[str]) -> np.ndarray:

    parsed: List[float] = []
    for item in values:
        try:
            parsed.append(float(item))
        except (TypeError, ValueError):
            raise ValueError(f"RCW CSV contains a nonnumeric value: {item!r}") from None
    return np.asarray(parsed, dtype=float)


def _series_type(raw_label: str, fallback_index: int, total_series: int) -> str:

    label = (raw_label or "").strip().lower()
    if "alpha" in label or "α" in label:
        return "alpha"
    if "beta" in label or "β" in label:
        return "beta"
    if not raw_label and total_series == 2:
        return "alpha" if fallback_index == 0 else "beta"
    return ""


def plot_csv(ax, csv_path: str) -> None:

    header, data = _read_csv(csv_path)


    if [label.strip().lower() for label in header] == ["col1", "col2", "col3"]:
        header = ["lambda", "alpha", "beta"]



    x_idx = _column_by_hint(header, ["lambda", "wavelength", "λ"])
    if x_idx < 0:
        x_idx = _column_by_hint(header, ["order", "index"])
    if x_idx < 0:
        x_idx = next((idx for idx, name in enumerate(header) if name.strip().lower() == "m"), 0)

    x_vals = _to_float([row[x_idx] if x_idx < len(row) else "" for row in data])
    if not np.all(np.isfinite(x_vals)):
        raise ValueError("RCW CSV contains a nonfinite scan coordinate.")


    hint_priority = [
        "alpha", "α", "beta", "β", "diff", "eff", "rs", "rp", "s11", "intensity", "power",
    ]
    y_indices: List[int] = []
    for hint in hint_priority:
        idx = _column_by_hint(header, [hint])
        if idx >= 0 and idx != x_idx and idx not in y_indices:
            y_indices.append(idx)
    if not y_indices:
        y_indices = [idx for idx in range(len(header)) if idx != x_idx]

    if not y_indices:
        raise ValueError("No value columns detected for plotting.")

    ax.cla()
    plotted_types: List[str] = []

    for series_pos, idx in enumerate(y_indices):
        series = _to_float([row[idx] if idx < len(row) else "" for row in data])
        mask = np.isfinite(x_vals) & np.isfinite(series)
        if not np.any(mask):
            continue

        raw_label = header[idx].strip()
        normalized_key = _series_type(raw_label, series_pos, len(y_indices))

        if normalized_key:
            label = _GREEK_LABELS.get(normalized_key, raw_label or f"col{idx + 1}")
        else:
            label = raw_label or f"col{idx + 1}"

        style = {"linewidth": 2.5}
        if np.count_nonzero(mask) == 1:
            style["marker"] = "o"
        if normalized_key and normalized_key in _GREEK_COLORS:
            style["color"] = _GREEK_COLORS[normalized_key]

            style["linestyle"] = "--" if normalized_key == "beta" else "-"

        ax.plot(x_vals, np.where(np.isfinite(series), series, np.nan), label=label, **style)

        if normalized_key:
            plotted_types.append(normalized_key)

    if not ax.lines:
        raise ValueError("Parsed CSV but found no plottable numeric data.")


    xlabel = header[x_idx].strip() or "Diffraction Order"
    normalized_xlabel = xlabel.lower()
    is_wavelength = any(hint in normalized_xlabel for hint in ("lambda", "wave", "λ"))

    ax.set_xlabel("Wavelength λ (μm)" if is_wavelength else xlabel)
    ax.set_ylabel("Normalized Stokes parameter" if plotted_types else "Value")

    finite_x = x_vals[np.isfinite(x_vals)]
    if is_wavelength and finite_x.size:
        span = float(finite_x.max() - float(finite_x.min()))
        step = _nice_step(span / 6.0) if span else 0.1
        ax.xaxis.set_major_locator(MultipleLocator(step))
        ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    else:
        ax.xaxis.set_major_locator(MultipleLocator(1))
        ax.xaxis.set_minor_locator(AutoMinorLocator(2))

    ax.grid(True, which="major", linewidth=0.8, alpha=0.6)
    ax.grid(True, which="minor", linewidth=0.4, alpha=0.35)

    if ax.lines:
        ax.legend(loc="best", frameon=False)

    if {"alpha", "beta"}.issubset(set(plotted_types)):
        x_title = r"Wavelength $\lambda$" if is_wavelength else xlabel
        ax.set_title(r"$\alpha$ and $\beta$ vs " + x_title, fontsize=16, fontstyle="italic")
    else:
        ax.set_title("")

    ax.figure.tight_layout()
