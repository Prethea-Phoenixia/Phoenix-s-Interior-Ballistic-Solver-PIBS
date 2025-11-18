__version__ = "0.5.5"


import logging

# import sys
from typing import Any

DESCRIPTION = "Description"
FONTNAME = "Sarasa Fixed SC"
FONTSIZE = 9
BOLDSIZE = 10
CONTEXT: dict[str, Any] = {
    "font.size": FONTSIZE,
    "axes.titlesize": FONTSIZE,
    "axes.labelsize": FONTSIZE,
    "axes.titlelocation": "right",
    "xtick.labelsize": FONTSIZE,
    "ytick.labelsize": FONTSIZE,
    "legend.fontsize": FONTSIZE,
    "figure.titlesize": BOLDSIZE,
    "lines.markersize": FONTSIZE / 4,
    "lines.linewidth": 1,
    "font.family": FONTNAME,
    "axes.labelweight": "bold",
    "xaxis.labellocation": "right",
    "yaxis.labellocation": "top",
}


THEMES = {
    "awlight": {"cmap": "afmhot_r", "is_light": True},
    "awdark": {"cmap": "afmhot", "is_light": False},
    "awwinxpblue": {"cmap": "afmhot_r", "is_light": True},
    "awclearlooks": {"cmap": "afmhot_r", "is_light": True},
}


root_logger = logging.getLogger(__name__)
root_logger.setLevel(logging.INFO)
log_formatter = logging.Formatter("%(asctime)s %(module)s.%(funcName)s: %(message)s")
# handler = logging.StreamHandler(sys.stderr)
# handler.setFormatter(formatter)
# root_logger.addHandler(handler)
