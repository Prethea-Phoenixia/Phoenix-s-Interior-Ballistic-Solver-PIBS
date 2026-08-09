__version__ = "0.5.5"


import logging
import sys

DESCRIPTION = "Description"
FONTNAME = "Sarasa Fixed SC"
FONTSIZE = 9
BOLDSIZE = 9


THEMES = {
    "awlight": {"cmap": "afmhot_r", "is_light": True},
    "awdark": {"cmap": "afmhot", "is_light": False},
}


root_logger = logging.getLogger(__name__)
root_logger.setLevel(logging.INFO)
log_formatter = logging.Formatter("%(asctime)s.%(msecs)03d - %(message)s", datefmt="%H:%M:%S")
handler = logging.StreamHandler(sys.stderr)
handler.setFormatter(log_formatter)
root_logger.addHandler(handler)
