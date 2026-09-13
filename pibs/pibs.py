import logging
import multiprocessing

from tkinter import Menu, Tk
from tkinter.font import Font

import matplotlib as mpl
from matplotlib import font_manager

from . import FONTNAME, FONTSIZE, __version__
from .ib_frame import InteriorBallisticsFrame
from .misc import (
    STRING,
    detect_darkmode_in_windows,
    get_windows_locale,
    loadfont,
    resolve_path,
    setup_windows_dpi,
    unloadfont,
)

logger = logging.getLogger(__name__)


class PIBS(Tk):
    def __init__(self, *args, loc: str, debug: bool, **kwargs):
        super().__init__(*args, **kwargs)
        # Only works on Windows:(
        if self._is_windows():
            super().iconbitmap(default=resolve_path("ui/logo.ico"))
        super().option_add("*tearOff", False)

        self.dpi: float = self.winfo_fpixels("1i")
        self.tk.call("tk", "scaling", 1.0 * self.dpi / 72.0)

        font = Font(family=FONTNAME, size=FONTSIZE)
        self.option_add("*Font", font)

        self.tk.call("lappend", "auto_path", resolve_path("ui/tksvg0.14"))
        self.tk.call("lappend", "auto_path", resolve_path("ui/awthemes-10.4.0"))

        self.title("PIBS v" + __version__)
        menubar = Menu(self)
        self.config(menu=menubar)

        self.bind("<Escape>", lambda event: self.state("normal"))
        self.bind("<F11>", lambda event: self.state("zoomed"))

        self.rowconfigure(0, weight=1)
        self.columnconfigure(0, weight=1)
        self.frame = InteriorBallisticsFrame(
            self,
            menubar,
            default_lang="English" if loc != "zh_CN" else "中文",
            localization_dict=STRING,
            font=font,
            os_dark=detect_darkmode_in_windows(),
            debug=debug,
        )
        self.frame.grid(row=0, column=0, sticky="nsew")

        self.minsize(self.winfo_width(), self.winfo_height())  # set minimum size
        self.is_fullscreen = False

        self.bind("<F4>", lambda *_: self.toggle_fullscreen())

    @staticmethod
    def _is_windows() -> bool:
        import platform

        return platform.system() == "Windows"

    def toggle_fullscreen(self):
        self.wm_attributes("-fullscreen", self.is_fullscreen)
        self.is_fullscreen = not self.is_fullscreen

    def quit(self):
        # explicitly unload the font at the end of program.
        unloadfont(resolve_path("ui/SarasaFixedSC-Regular.ttf"), True, True)
        super().quit()


def main(loc: str = "", debug: bool = False):
    multiprocessing.freeze_support()
    logger.info("Initializing")

    setup_windows_dpi()
    loc = loc or get_windows_locale()

    loadfont(resolve_path("ui/SarasaFixedSC-Regular.ttf"), True, True)
    mpl.font_manager.fontManager.addfont(resolve_path("ui/SarasaFixedSC-Regular.ttf"))

    pibs = PIBS(loc=loc, debug=debug)
    pibs.mainloop()


if __name__ == "__main__":
    main()
