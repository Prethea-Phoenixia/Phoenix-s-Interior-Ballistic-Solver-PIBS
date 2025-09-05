import json
import os
import sys
from ctypes import byref, create_string_buffer, create_unicode_buffer, windll
from math import floor, log, log10

_prefix = {
    "y": 1e-24,  # yocto
    "z": 1e-21,  # zepto
    "a": 1e-18,  # atto
    "f": 1e-15,  # femto
    "p": 1e-12,  # pico
    "n": 1e-9,  # nano
    "μ": 1e-6,  # micro
    "m": 1e-3,  # mili
    # "c": 1e-2,  # centi
    # "d": 1e-1,  # deci
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


# noinspection PyTypeChecker
def resolvepath(path):
    if getattr(sys, "frozen", False):
        # If the 'frozen' flag is set, we are in bundled-app mode!
        resolved_path = os.path.abspath(os.path.join(sys._MEIPASS, path))
    else:
        # Normal development mode. Use os.getcwd() or __file__ as appropriate in your case...
        resolved_path = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), path))

    return resolved_path


FR_PRIVATE = 0x10
FR_NOT_ENUM = 0x20


def loadfont(fontpath, private=True, enumerable=False):
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
    if isinstance(fontpath, bytes):
        pathbuf = create_string_buffer(fontpath)
        add_font_resource_ex = windll.gdi32.AddFontResourceExA

    elif isinstance(fontpath, str):
        pathbuf = create_unicode_buffer(fontpath)
        add_font_resource_ex = windll.gdi32.AddFontResourceExW
    else:
        raise TypeError("fontpath must be of type str or unicode")

    flags = (FR_PRIVATE if private else 0) | (FR_NOT_ENUM if not enumerable else 0)
    num_fonts_added = add_font_resource_ex(byref(pathbuf), flags, 0)
    return bool(num_fonts_added)


def unloadfont(fontpath, private=True, enumerable=False):
    """
    Unloads the fonts in the specified file.

    see http://msdn2.microsoft.com/en-us/library/ms533925.aspx
    """
    if isinstance(fontpath, bytes):
        pathbuf = create_string_buffer(fontpath)
        remove_font_resource_ex = windll.gdi32.RemoveFontResourceExA
    elif isinstance(fontpath, str):
        pathbuf = create_unicode_buffer(fontpath)
        remove_font_resource_ex = windll.gdi32.RemoveFontResourceExW
    else:
        raise TypeError("fontpath must be a str or unicode")

    flags = (FR_PRIVATE if private else 0) | (FR_NOT_ENUM if not enumerable else 0)
    return bool(remove_font_resource_ex(byref(pathbuf), flags, 0))


def to_si(v, dec=4, unit=None, use_sn=False):
    if v is None:
        return "N/A"
    elif isinstance(v, int) or isinstance(v, float):
        pass
    else:
        raise ValueError(f"Cannot convert type of {type(v):} to SI notation")

    if v >= 0:
        positive = True
    else:
        positive = False
    v = abs(v)
    for prefix, magnitude, nextMagnitude in zip(
        _prefix.keys(),
        tuple(_prefix.values())[:-1],
        tuple(_prefix.values())[1:],
    ):
        if 1 <= (v / magnitude) < (nextMagnitude / magnitude):
            # allows handling of non-uniformly log10 spaced prefixes
            vstr = "{:#.{:}g}".format(v / magnitude, dec)
            return (
                (" " if positive else "-")
                + vstr
                + " " * (dec + 1 - len(vstr) + vstr.find("."))
                + (
                    ("E{:<+3d}".format(round(log(magnitude, 10))) if magnitude != 1 else "    ")  # 4 SPACES!
                    if use_sn
                    else prefix
                )
                + (unit if unit is not None else "")
            )
    if v == 0:
        return (
            (" " if positive else "-")
            + "{:#.{:}g}".format(v, dec)
            + ("     " if use_sn else "  ")
            + (unit if unit is not None else "")
        )
    else:  # return a result in SI as a last resort
        closest = log(v, 10) // 3
        magnitude = 10 ** (closest * 3)
        vstr = "{:#.{:}g}".format(v / magnitude, dec)
        return (
            (" " if positive else "-")
            + vstr
            + " " * (dec + 1 - len(vstr) + vstr.find("."))
            + ("E{:<+3d}".format(round(log(magnitude, 10))))
            + (unit if unit is not None else "")
        )


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


def dot_aligned(matrix, units, use_sn, strip_ws=True):
    transposed = []

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


def round_sig(x, n=4):
    return round(x, (n - 1) - int(floor(log10(abs(x)))))


def format_mass(m, n=4):
    if m:
        if m < 1e-3:
            return "{:.{:}g} mg".format(m * 1e6, n)
        elif m < 1:
            return "{:.{:}g} g".format(m * 1e3, n)
        elif m < 1000:
            return "{:.{:}g} kg".format(m, n)
        elif m < 1e6:
            return "{:.{:}g} t".format(m * 1e-3, n)
        elif m < 1e9:
            return "{:.{:}g} kt".format(m * 1e-6, n)
        elif m < 1e12:
            return "{:.{:}g} Mt".format(m * 1e-9, n)

    return "N/A"


with open(resolvepath("ui/localization.json"), encoding="utf-8") as file:
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
