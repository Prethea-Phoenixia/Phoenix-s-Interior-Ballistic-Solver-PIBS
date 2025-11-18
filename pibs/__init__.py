__version__ = "0.5.5"


import logging
import sys
from typing import Any

DESCRIPTION = "Description"
USE_LF = "USE_LF"
USE_CV = "USE_CV"
USE_LD = "USE_LD"

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


root_logger = logging.getLogger()
root_logger.setLevel(logging.INFO)
handler = logging.StreamHandler(sys.stderr)
formatter = logging.Formatter("%(asctime)s - %(name)s - %(message)s")
handler.setFormatter(formatter)
root_logger.addHandler(handler)


class TextHandler(logging.Handler):
    # This class allows you to log to a Tkinter Text or ScrolledText widget
    # Adapted from Moshe Kaplan: https://gist.github.com/moshekaplan/c425f861de7bbf28ef06

    def __init__(self, text):
        # run the regular Handler __init__
        logging.Handler.__init__(self)
        self.setFormatter(formatter)
        self.setLevel(logging.INFO)
        # Store a reference to the Text it will log to
        self.text = text
        self.text.config(state="disabled")

    def emit(self, record):
        msg = self.format(record)

        def append():
            self.text.configure(state="normal")

            # use record name as tags
            tags = record.name.split(".")
            self.text.insert("end", msg.strip("\n") + "\n", [record.levelno, *tags])
            self.text.configure(state="disabled")

            # Autoscroll to the bottom
            self.text.yview("end")

        # This is necessary because we can't modify the Text from other threads
        # self.text.after(1000, append)
        self.text.after_idle(append)
