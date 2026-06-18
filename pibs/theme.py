from __future__ import annotations

from tkinter import ttk, StringVar
from . import THEMES, FONTNAME, BOLDSIZE, FONTSIZE, DESCRIPTION


class ThemedMixin(object):
    def __init__(self, master, *args, dpi: int, os_dark: bool = True, **kwargs):
        super().__init__(master, *args, **kwargs)

        if isinstance(master, ThemedMixin):
            self.theme_name_var = master.theme_name_var
        else:
            self.theme_name_var = StringVar(value="awdark" if os_dark else "awlight")
        self.dpi = dpi
        self.force_update_on_theme_widget = []

        self.context = {
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

    def use_theme(self):
        style = ttk.Style(self)
        style.theme_use(self.theme_name_var.get())

        """ensure that the treeview rows are roughly the same height
        regardless of dpi. on Windows, default is Segoe UI at 9 points
        so the default row height should be around 12"""

        style.configure("Treeview", rowheight=round(12 * (FONTSIZE / 9) * self.dpi / 72.0))
        style.configure("Treeview.Heading", font=(FONTNAME, FONTSIZE))
        style.configure("TButton", font=(FONTNAME, BOLDSIZE, "bold"))
        style.configure("TLabelframe.Label", font=(FONTNAME, BOLDSIZE, "bold"))
        style.configure("TCheckbutton", font=(FONTNAME, FONTSIZE))
        style.configure("TNotebook.Tab", font=(FONTNAME, BOLDSIZE, "bold"))

        bgc = str(style.lookup("TEntry", "background"))
        fgc = str(style.lookup("TEntry", "foreground"))
        fbgc = str(style.lookup("TEntry", "fieldbackground")) or bgc

        # some widgets also needs to be manually updated
        for widget in self.force_update_on_theme_widget:
            widget.config(background=fbgc, foreground=fgc, insertbackground=fgc)

        self.context.update(
            {
                "figure.facecolor": bgc,
                "figure.edgecolor": fgc,
                "axes.edgecolor": fgc,
                "axes.labelcolor": fgc,
                "axes.facecolor": fbgc,
                "text.color": fgc,
                "xtick.color": fgc,
                "ytick.color": fgc,
            }
        )
