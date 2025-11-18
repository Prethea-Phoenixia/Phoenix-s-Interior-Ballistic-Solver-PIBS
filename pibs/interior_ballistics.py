import locale

import multiprocessing
import platform

if platform.system() == "Windows":
    from ctypes import windll

from tkinter import Menu, Tk
from tkinter.font import Font

import matplotlib as mpl

from matplotlib import font_manager

from . import __version__

from .misc import STRING, detect_darkmode_in_windows, loadfont, resolve_path, unloadfont

from . import FONTSIZE, FONTNAME
from .interior_ballistics_frame import InteriorBallisticsFrame
from . import root_logger


def grid_configure_recursive(widget, **kwargs):
    stack = list(widget.winfo_children())
    while stack:
        descendent = stack.pop()
        stack.extend(descendent.winfo_children())
        descendent.grid_configure(**kwargs)


class PIBS(Tk):
    def __init__(self, *args, loc: str, debug: bool, **kwargs):
        super().__init__(*args, **kwargs)
        # Only works on Windows:(
        if platform.system() == "Windows":
            super().iconbitmap(default=resolve_path("ui/logo.ico"))
        super().option_add("*tearOff", False)

        font = Font(family=FONTNAME, size=FONTSIZE)
        self.option_add("*Font", font)

        dpi = self.winfo_fpixels("1i")
        # Tk was originally developed for a dpi of 72
        self.tk.call("tk", "scaling", 1.0 * dpi / 72.0)
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
            dpi,
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

    def toggle_fullscreen(self):
        self.wm_attributes("-fullscreen", self.is_fullscreen)
        self.is_fullscreen = not self.is_fullscreen

    def quit(self):
        # explicitly unload the font at the end of program.
        unloadfont(resolve_path("ui/SarasaFixedSC-Regular.ttf"), True, True)
        super().quit()


def main(loc: str = None, debug: bool = False):
    multiprocessing.freeze_support()
    root_logger.info("Initializing")

    # this tells windows that our program will handle scaling ourselves
    if platform.system() == "Windows":
        win_release = platform.release()
        if win_release in ("8", "10", "11"):
            # noinspection PyUnresolvedReferences
            windll.shcore.SetProcessDpiAwareness(1)
        elif win_release in ("7", "Vista"):
            # noinspection PyUnresolvedReferences
            windll.user32.SetProcessDPIAware()

        # noinspection PyUnresolvedReferences
        loc = loc if loc else locale.windows_locale[windll.kernel32.GetUserDefaultUILanguage()]

    loadfont(resolve_path("ui/SarasaFixedSC-Regular.ttf"), True, True)
    mpl.font_manager.fontManager.addfont(resolve_path("ui/SarasaFixedSC-Regular.ttf"))

    pibs = PIBS(loc=loc, debug=debug)
    pibs.mainloop()


if __name__ == "__main__":
    main()
