from __future__ import annotations

import json
import logging
import os
import platform
import shutil
import sys
from ctypes import byref, create_string_buffer, create_unicode_buffer
from pathlib import Path

logger = logging.getLogger(__name__)

if platform.system() == "Windows":
    from ctypes import windll

from math import floor, log, log10, isinf, isnan

_superscript_map = {
    "0": "⁰",
    "1": "¹",
    "2": "²",
    "3": "³",
    "4": "⁴",
    "5": "⁵",
    "6": "⁶",
    "7": "⁷",
    "8": "⁸",
    "9": "⁹",
    "a": "ᵃ",
    "b": "ᵇ",
    "c": "ᶜ",
    "d": "ᵈ",
    "e": "ᵉ",
    "f": "ᶠ",
    "g": "ᵍ",
    "h": "ʰ",
    "i": "ᶦ",
    "j": "ʲ",
    "k": "ᵏ",
    "l": "ˡ",
    "m": "ᵐ",
    "n": "ⁿ",
    "o": "ᵒ",
    "p": "ᵖ",
    "q": "۹",
    "r": "ʳ",
    "s": "ˢ",
    "t": "ᵗ",
    "u": "ᵘ",
    "v": "ᵛ",
    "w": "ʷ",
    "x": "ˣ",
    "y": "ʸ",
    "z": "ᶻ",
    "A": "ᴬ",
    "B": "ᴮ",
    "C": "ᶜ",
    "D": "ᴰ",
    "E": "ᴱ",
    "F": "ᶠ",
    "G": "ᴳ",
    "H": "ᴴ",
    "I": "ᴵ",
    "J": "ᴶ",
    "K": "ᴷ",
    "L": "ᴸ",
    "M": "ᴹ",
    "N": "ᴺ",
    "O": "ᴼ",
    "P": "ᴾ",
    "Q": "Q",
    "R": "ᴿ",
    "S": "ˢ",
    "T": "ᵀ",
    "U": "ᵁ",
    "V": "ⱽ",
    "W": "ᵂ",
    "X": "ˣ",
    "Y": "ʸ",
    "Z": "ᶻ",
    "+": "⁺",
    "-": "⁻",
    "=": "⁼",
    "(": "⁽",
    ")": "⁾",
}

_subscript_map = {
    "0": "₀",
    "1": "₁",
    "2": "₂",
    "3": "₃",
    "4": "₄",
    "5": "₅",
    "6": "₆",
    "7": "₇",
    "8": "₈",
    "9": "₉",
    "a": "ₐ",
    "b": "♭",
    "c": "꜀",
    "d": "ᑯ",
    "e": "ₑ",
    "f": "բ",
    "g": "₉",
    "h": "ₕ",
    "i": "ᵢ",
    "j": "ⱼ",
    "k": "ₖ",
    "l": "ₗ",
    "m": "ₘ",
    "n": "ₙ",
    "o": "ₒ",
    "p": "ₚ",
    "q": "૧",
    "r": "ᵣ",
    "s": "ₛ",
    "t": "ₜ",
    "u": "ᵤ",
    "v": "ᵥ",
    "w": "w",
    "x": "ₓ",
    "y": "ᵧ",
    "z": "₂",
    "A": "ₐ",
    "B": "₈",
    "C": "C",
    "D": "D",
    "E": "ₑ",
    "F": "բ",
    "G": "G",
    "H": "ₕ",
    "I": "ᵢ",
    "J": "ⱼ",
    "K": "ₖ",
    "L": "ₗ",
    "M": "ₘ",
    "N": "ₙ",
    "O": "ₒ",
    "P": "ₚ",
    "Q": "Q",
    "R": "ᵣ",
    "S": "ₛ",
    "T": "ₜ",
    "U": "ᵤ",
    "V": "ᵥ",
    "W": "w",
    "X": "ₓ",
    "Y": "ᵧ",
    "Z": "Z",
    "+": "₊",
    "-": "₋",
    "=": "₌",
    "(": "₍",
    ")": "₎",
}

sup_trans = str.maketrans("".join(_superscript_map.keys()), "".join(_superscript_map.values()))

sub_trans = str.maketrans("".join(_subscript_map.keys()), "".join(_subscript_map.values()))

_prefix = {
    "y": 1e-24,  # yocto
    "z": 1e-21,  # zepto
    "a": 1e-18,  # atto
    "f": 1e-15,  # femto
    "p": 1e-12,  # pico
    "n": 1e-9,  # nano
    "μ": 1e-6,  # micro
    "m": 1e-3,  # mili
    "c": 1e-2,  # centi
    "d": 1e-1,  # deci
    " ": 1,  # unit
    # "da": 1e1, # deca
    # "h": 1e2,  # hecto
    "k": 1e3,  # kilo
    "M": 1e6,  # mega
    "G": 1e9,  # giga
    "T": 1e12,  # tera
    "P": 1e15,  # peta
    "E": 1e18,  # exa
    "Z": 1e21,  # zetta
    "Y": 1e24,  # yotta
}


def get_font_dir() -> Path | None:
    if platform.system() == "Windows":
        return None
    else:
        font_dir = Path.home() / ".local" / "share" / "fonts"
        font_dir.mkdir(parents=True, exist_ok=True)
        return font_dir


def resolve_path(path: str) -> str:
    if getattr(sys, "frozen", False):
        # If the 'frozen' flag is set, we are in bundled-app mode!
        resolved_path = os.path.abspath(os.path.join(sys._MEIPASS, path))
    else:
        # Normal development mode. Use os.getcwd() or __file__ as appropriate in your case...
        resolved_path = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), path))

    return resolved_path


FR_PRIVATE = 0x10
FR_NOT_ENUM = 0x20


def loadfont(fontpath, private: bool = True, enumerable: bool = False) -> bool:
    """
    Makes fonts located in file `fontpath` available to the font system.
        `private`     if True, other processes cannot see this font, and this
                      font will be unloaded when the process dies
        `enumerable`  if True, this font will appear when enumerating fonts

    See https://msdn.microsoft.com/en-us/library/dd183327(VS.85).aspx

    This function was taken from
    https://github.com/ifwe/digsby/blob/f5fe00244744aa131e07f09348d10563f3d8fa99/digsby/src/gui/native/win/winfonts.py#L15
    This function is written for Python 2.x. For 3.x, you
    have to convert the isinstance checks to bytes and str
    """
    if platform.system() == "Windows":
        if isinstance(fontpath, bytes):
            pathbuf = create_string_buffer(fontpath)
            # noinspection PyUnresolvedReferences
            add_font_resource_ex = windll.gdi32.AddFontResourceExA
        elif isinstance(fontpath, str):
            pathbuf = create_unicode_buffer(fontpath)
            # noinspection PyUnresolvedReferences
            add_font_resource_ex = windll.gdi32.AddFontResourceExW
        else:
            raise TypeError("fontpath must be of type str or unicode")

        flags = (FR_PRIVATE if private else 0) | (FR_NOT_ENUM if not enumerable else 0)
        num_fonts_added = add_font_resource_ex(byref(pathbuf), flags, 0)
        logger.info(f"Windows: Loaded font {fontpath} (private={private}, enumerable={enumerable})")
        return bool(num_fonts_added)
    else:
        font_dir = get_font_dir()
        if font_dir is None:
            return False
        shutil.copy2(fontpath, font_dir)
        logger.info(f"Linux: Copied font {fontpath} to {font_dir}")
        return True


def unloadfont(fontpath, private: bool = True, enumerable: bool = False) -> bool:
    """
    Unloads the fonts in the specified file.

    see http://msdn2.microsoft.com/en-us/library/ms533925.aspx
    """
    if platform.system() == "Windows":
        if isinstance(fontpath, bytes):
            pathbuf = create_string_buffer(fontpath)
            # noinspection PyUnresolvedReferences
            remove_font_resource_ex = windll.gdi32.RemoveFontResourceExA
        elif isinstance(fontpath, str):
            pathbuf = create_unicode_buffer(fontpath)
            # noinspection PyUnresolvedReferences
            remove_font_resource_ex = windll.gdi32.RemoveFontResourceExW
        else:
            raise TypeError("fontpath must be a str or unicode")

        flags = (FR_PRIVATE if private else 0) | (FR_NOT_ENUM if not enumerable else 0)
        logger.info(f"Windows: Unloaded font {fontpath} (private={private}, enumerable={enumerable})")
        return bool(remove_font_resource_ex(byref(pathbuf), flags, 0))
    else:
        font_dir = get_font_dir()
        if font_dir is None:
            return False
        dest_path = font_dir / Path(fontpath).name
        if dest_path.exists():
            dest_path.unlink()
        logger.info(f"Linux: Removed font {fontpath} from {font_dir}")
        return True


def to_si(v: float, dec: int = 4, unit: str = "", unit_dim: int = 1, use_sn: bool = False) -> str:
    if v is None:
        return "N/A"
    elif isinstance(v, int) or isinstance(v, float):
        if isinf(v):
            return "INF"
        if isnan(v):
            return "NAN"
        if v == 0:
            return " " + "{:#.{:}g}".format(v, dec) + ("     " if use_sn else "  ") + (unit if unit is not None else "")
        else:
            if v > 0:
                positive = True
            else:
                positive = False

            v = abs(v)
            for prefix, magnitude, next_magnitude in zip(
                _prefix.keys(),
                tuple(_prefix.values())[:-1],
                tuple(_prefix.values())[1:],
            ):
                if 1 <= (v / (magnitude**unit_dim)) < (next_magnitude / (magnitude**unit_dim)):
                    # allows handling of non-uniformly log10 spaced prefixes
                    vstr = "{:#.{:}g}".format(v / magnitude**unit_dim, dec)
                    return (
                        (" " if positive else "-")
                        + vstr
                        + " " * (dec + 1 - len(vstr) + vstr.find("."))
                        + (("×⏨" + f"{int(log(magnitude, 10)): <3d}".translate(sup_trans)) if use_sn else prefix)
                        + (unit if unit is not None else "")
                    )

        raise ValueError(f"Cannot convert {v} to SI notation")

    elif isinstance(v, str):
        return v
    else:
        raise ValueError(f"Cannot convert type of {type(v):} to SI notation")


def validate_nn(inp):
    """
    validate an input if it results in:
    - result >=0
    - result is empty
    in the latter case, the empty field will be handled by tracing
    change in variable.
    """
    if (
        inp == ""
        or inp == "."
        or (inp.count("e") == 1 and inp[-1] == "e")  # scientific input
        or (inp.count("e") == 1 and inp[-2:] == "e-")  # scientific input with negative exponent
    ):
        return True
    try:
        if float(inp) >= 0:
            return True
        else:
            return False
    except ValueError:
        return False


def validate_pi(inp):  # validate an input such that the result is a positive integer
    if inp == "":
        return True  # we will catch this by filling the default value
    try:
        return float(inp).is_integer() and float(inp) > 0 and "." not in inp
    except ValueError:
        return False


def validate_ce(inp: float) -> float:
    return validate_range(inp, low=0.0, high=100.0)


def validate_range(inp: float, low: float, high: float) -> float:  # validate a range
    high, low = max((high, low)), min((high, low))
    if inp == "":
        return True
    try:
        return high >= float(inp) >= low
    except ValueError:
        return False


def validate_flt(inp):  # validate an input such that the result is a float.
    if (
        inp == ""
        or inp == "."
        or inp == "-"
        or (inp.count("e") == 1 and inp[-1] == "e")  # scientific input
        or (inp.count("e") == 1 and inp[-2:] == "e-")  # scientific input with negative exponent
    ):
        return True
    try:
        float(inp)
        return True
    except ValueError:
        return False


def format_float_input(event, var):
    # v = event.widget.get()
    v = var.get()
    if v == "" or v == ".":
        var.set(0.0)
    else:
        var.set(float(v))


def format_int_input(event, var):
    # v = event.widget.get()
    v = var.get()
    if v == "":
        var.set(0)
    else:
        var.set(int(v))


def dot_aligned(matrix, units: tuple[str, ...], use_sn: tuple[bool, ...], strip_ws: bool = True):
    transposed = []  # TODO: properly annotate this thing.

    for seq, unit, isSN in zip(zip(*matrix), units, use_sn):
        snums = []
        for n in seq:
            if isinstance(n, int) or isinstance(n, float) or n is None:
                vstr = to_si(n, unit=unit, use_sn=isSN)
                snums.append(vstr.strip() if strip_ws else vstr)
            elif isinstance(n, str):
                snums.append(n)
            else:
                raise ValueError("Unknown type encountered in dot_aligned")
        dots = [s.find(".") for s in snums]
        m = max(dots)
        transposed.append(tuple(" " * (m - d) + s for s, d in zip(snums, dots)))

    return tuple(zip(*transposed))


def round_sig(x: float, n: int = 4):
    return round(x, (n - 1) - int(floor(log10(abs(x)))))


def format_mass(m: float, n: int = 4) -> str:
    if m:
        if m < 1e-6:
            return "{:.{:}g} μg".format(m * 1e9, n)
        elif m < 1e-3:
            return "{:.{:}g} mg".format(m * 1e6, n)
        elif m < 1:
            return "{:.{:}g}  g".format(m * 1e3, n)
        elif m < 1000:
            return "{:.{:}g} kg".format(m, n)
        elif m < 1e6:
            return "{:.{:}g}  t".format(m * 1e-3, n)
        elif m < 1e9:
            return "{:.{:}g} kt".format(m * 1e-6, n)
        elif m < 1e12:
            return "{:.{:}g} Mt".format(m * 1e-9, n)

    return "N/A"


def format_temp(t: float) -> str:
    if t:
        return f"{t:.0f} K "
    return "N/A"


with open(resolve_path("ui/localization.json"), encoding="utf-8") as file:
    STRING = json.load(file)


def filenameize(string: str) -> str:
    return string.replace(" ", "_")


def detect_darkmode_in_windows() -> bool:
    try:
        import winreg
    except ImportError:
        return False
    registry = winreg.ConnectRegistry(None, winreg.HKEY_CURRENT_USER)
    reg_keypath = r"SOFTWARE\Microsoft\Windows\CurrentVersion\Themes\Personalize"
    try:
        reg_key = winreg.OpenKey(registry, reg_keypath)
    except FileNotFoundError:
        return False

    for i in range(1024):
        try:
            value_name, value, _ = winreg.EnumValue(reg_key, i)
            if value_name == "AppsUseLightTheme":
                return value == 0
        except OSError:
            break
    return False


if __name__ == "__main__":
    print(to_si(1e-4).strip())
    from math import pi

    print(round_sig(pi))
