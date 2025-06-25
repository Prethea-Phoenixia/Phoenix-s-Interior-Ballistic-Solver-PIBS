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


def resolvepath(path):
    if getattr(sys, "frozen", False):
        # If the 'frozen' flag is set, we are in bundled-app mode!
        resolved_path = os.path.abspath(os.path.join(sys._MEIPASS, path))
    else:
        # Normal development mode. Use os.getcwd() or __file__ as appropriate in your case...
        resolved_path = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), path))

    return resolved_path


def loadfont(fontpath, private=True, enumerable=False):
    # Makes fonts located in file `fontpath` available to the font system.
    # `private`     if True, other processes cannot see this font, and this
    #               font will be unloaded when the process dies
    # `enumerable`  if True, this font will appear when enumerating fonts
    #
    # See https://msdn.microsoft.com/en-us/library/dd183327(VS.85).aspx

    # This function was taken from
    # https://github.com/ifwe/digsby/blob/f5fe00244744aa131e07f09348d10563f3d8fa99/digsby/src/gui/native/win/winfonts.py#L15
    # This function is written for Python 2.x. For 3.x, you
    # have to convert the isinstance checks to bytes and str

    FR_PRIVATE = 0x10
    FR_NOT_ENUM = 0x20
    if isinstance(fontpath, bytes):
        pathbuf = create_string_buffer(fontpath)
        AddFontResourceEx = windll.gdi32.AddFontResourceExA

    elif isinstance(fontpath, str):
        pathbuf = create_unicode_buffer(fontpath)
        AddFontResourceEx = windll.gdi32.AddFontResourceExW
    else:
        raise TypeError("fontpath must be of type str or unicode")

    flags = (FR_PRIVATE if private else 0) | (FR_NOT_ENUM if not enumerable else 0)
    numFontsAdded = AddFontResourceEx(byref(pathbuf), flags, 0)
    return bool(numFontsAdded)


def toSI(v, dec=4, unit=None, useSN=False):
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
                    if useSN
                    else prefix
                )
                + (unit if unit is not None else "")
            )
    if v == 0:
        return (
            (" " if positive else "-")
            + "{:#.{:}g}".format(v, dec)
            + ("     " if useSN else "  ")
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


def validateNN(inp):
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


def validatePI(inp):  # validate an input such that the result is a positive integer
    if inp == "":
        return True  # we will catch this by filling the default value
    try:
        return float(inp).is_integer() and float(inp) > 0 and "." not in inp
    except ValueError:
        return False


def validateFLT(inp):  # validate an input such that the result is a float.
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


def formatFloatInput(event, var):
    # v = event.widget.get()
    v = var.get()
    if v == "" or v == ".":
        var.set(event.widget.default)
    else:
        var.set(float(v))


def formatIntInput(event, var):
    # v = event.widget.get()
    v = var.get()
    if v == "":
        var.set(event.widget.default)
    else:
        var.set(int(v))


def dot_aligned(matrix, units, useSN, stripWS=True):
    transposed = []

    for seq, unit, isSN in zip(zip(*matrix), units, useSN):
        snums = []
        for n in seq:
            if isinstance(n, int) or isinstance(n, float) or n is None:
                vstr = toSI(n, unit=unit, useSN=isSN)
                snums.append(vstr.strip() if stripWS else vstr)
            elif isinstance(n, str):
                snums.append(n)
            else:
                raise ValueError("Unknown type encountered in dot_aligned")
        dots = [s.find(".") for s in snums]
        m = max(dots)
        transposed.append(tuple(" " * (m - d) + s for s, d in zip(snums, dots)))

    return tuple(zip(*transposed))


def roundSig(x, n=4):
    return round(x, (n - 1) - int(floor(log10(abs(x)))))


def formatMass(m, n=4):
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


with open(resolvepath("ui/localization.json")) as file:
    STRING = json.load(file)


if __name__ == "__main__":
    print(toSI(1e-4).strip())
    from math import pi

    print(roundSig(pi))
