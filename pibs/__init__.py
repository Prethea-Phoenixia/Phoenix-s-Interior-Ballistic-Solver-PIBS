__version__ = "0.5.6"


import logging
import sys
from enum import Enum

DESCRIPTION = "Description"
FONTNAME = "Sarasa Fixed SC"
FONTSIZE = 9
BOLDSIZE = 10


class FigureType(str, Enum):
    MAIN = "main"
    AUX = "aux"
    GEOM = "geom"
    GUIDE = "guide"

    def __str__(self) -> str:
        return self.value


class Theme(str, Enum):
    AWDARK = "awdark"
    AWLIGHT = "awlight"

    def __str__(self) -> str:
        return self.value


THEMES = {
    Theme.AWLIGHT: {"cmap": "afmhot_r", "is_light": True},
    Theme.AWDARK: {"cmap": "afmhot", "is_light": False},
}


root_logger = logging.getLogger(__name__)
root_logger.setLevel(logging.INFO)
root_logger.propagate = False
log_formatter = logging.Formatter("%(asctime)s - %(message)s", datefmt="%H:%M:%S")

handler = logging.StreamHandler(sys.stderr)
handler.setFormatter(log_formatter)
root_logger.addHandler(handler)
