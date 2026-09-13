__version__ = "0.5.5"


import logging
import sys
from enum import Enum

DESCRIPTION = "Description"
FONTNAME = "Sarasa Fixed SC"
FONTSIZE = 9
BOLDSIZE = 9


class FigureType(str, Enum):
    MAIN = "main"
    AUX = "aux"
    GEOM = "geom"
    GUIDE = "guide"


class Theme(str, Enum):
    AWDARK = "awdark"
    AWLIGHT = "awlight"


THEMES = {
    Theme.AWLIGHT: {"cmap": "afmhot_r", "is_light": True},
    Theme.AWDARK: {"cmap": "afmhot", "is_light": False},
}


root_logger = logging.getLogger(__name__)
root_logger.setLevel(logging.INFO)
log_formatter = logging.Formatter("%(asctime)s - [%(module)10s] %(message)s", datefmt="%H:%M:%S")
handler = logging.StreamHandler(sys.stderr)
handler.setFormatter(log_formatter)
root_logger.addHandler(handler)
