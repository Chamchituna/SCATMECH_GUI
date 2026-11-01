import math
import os
import re
import shutil
from pathlib import Path

from scatmech_paths import get_scatmech_bin, native_path

_NUMBER = re.compile(r"[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?")


def single_line(value, label, *, allow_empty=False):
    if value is None or isinstance(value, (bool, dict, list, tuple)):
        raise ValueError(f"{label} must be text or a number.")
    text = str(value)
    if any(char in text for char in ("\n", "\r", "\x00")):
        raise ValueError(f"{label} must fit on one line.")
    text = text.strip()
    if not text and not allow_empty:
        raise ValueError(f"{label} is required.")
    return text


def finite_float(value, label, *, minimum=None, maximum=None, strict_minimum=False):
    text = single_line(value, label)
    if _NUMBER.fullmatch(text) is None:
        raise ValueError(f"{label} must be a finite number (decimal or scientific notation).")
    try:
        number = float(text)
    except (ValueError, OverflowError):
        raise ValueError(f"{label} must be a finite number.") from None
    if not math.isfinite(number):
        raise ValueError(f"{label} must be a finite number.")
    if minimum is not None and (number < minimum or (strict_minimum and number == minimum)):
        comparison = "greater than" if strict_minimum else "at least"
        raise ValueError(f"{label} must be {comparison} {minimum}.")
    if maximum is not None and number > maximum:
        raise ValueError(f"{label} must be at most {maximum}.")
    return number


def integer(value, label, *, minimum=None, maximum=None):
    lower = max(-2147483648, minimum) if minimum is not None else -2147483648
    upper = min(2147483647, maximum) if maximum is not None else 2147483647
    number = finite_float(value, label, minimum=lower, maximum=upper)
    if not number.is_integer():
        raise ValueError(f"{label} must be a whole number.")
    return int(number)


def optical_constant(value, label):

    text = single_line(value, label)
    if any(char in text for char in "(),"):
        pair_text = text[1:-1] if text.startswith("(") and text.endswith(")") else text
        pair = re.fullmatch(r"\s*([^(),]+)\s*,\s*([^(),]+)\s*", pair_text)
        if pair is None:
            raise ValueError(f"{label} must be an (n,k) pair, a real index, or a dielectric filename.")
        n = finite_float(pair.group(1), f"{label} n")
        k = finite_float(pair.group(2), f"{label} k")
        return f"({n:.15g},{k:.15g})"
    try:
        float(text)
    except ValueError:
        return solver_file(text, label)
    n = finite_float(text, label)
    return f"({n:.15g},0)"


def solver_file(value, label):

    text = single_line(value, label)
    roots = [Path.cwd(), get_scatmech_bin()]
    roots.extend(Path(part).expanduser() for part in os.environ.get("SCATMECH", "").split(os.pathsep) if part)
    for solver in ("brdfprog", "mieprog", "reflectprog", "rcwprog"):
        exe = shutil.which(solver)
        if exe:
            roots.append(Path(exe).parent)
    filename = Path(text).expanduser()
    for candidate in ([filename] if filename.is_absolute() else [root / filename for root in roots]):
        if candidate.is_file():
            resolved = str(native_path(candidate))
            if any(char.isspace() for char in resolved):
                raise ValueError(f"{label}: SCATMECH filenames cannot contain spaces; move the file to a path without spaces.")
            return resolved
    raise ValueError(f"{label}: file not found: {text}")
