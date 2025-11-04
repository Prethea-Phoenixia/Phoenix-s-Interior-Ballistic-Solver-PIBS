import csv
import json
import locale
import logging
import multiprocessing
import platform
import sys
import traceback
if platform.system() == "Windows":
    from ctypes import windll
from itertools import repeat
from logging.handlers import QueueHandler, QueueListener
from math import ceil, inf, log10
from multiprocessing import Process, Queue
from queue import Empty
from tkinter import (
    IntVar,
    Menu,
    StringVar,
    Text,
    Tk,
    filedialog,
    messagebox,
    ttk,
)
from tkinter.font import Font
from tkinter.ttk import Frame
from typing import Literal

import matplotlib as mpl
from labellines import labelLines
from matplotlib import font_manager
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

from pibs.guidegraph import guide_graph

from . import __version__
from .ballistics import (
    COMPUTE,
    CONVENTIONAL,
    DOMAIN_LEN,
    DOMAIN_TIME,
    POINT_BURNOUT,
    POINT_EXIT,
    POINT_FRACTURE,
    POINT_PEAK_AVG,
    POINT_PEAK_BREECH,
    POINT_PEAK_SHOT,
    POINT_PEAK_STAG,
    POINT_START,
    RECOILLESS,
    SOL_LAGRANGE,
    SOL_MAMONTOV,
    SOL_PIDDUCK,
    Composition,
    Constrained,
    ConstrainedRecoilless,
    Geometry,
    Gun,
    Material,
    Propellant,
    Recoilless,
    SimpleGeometry,
)
from .localized_widget import LocalizedFrame
from .misc import (
    STRING,
    detect_darkmode_in_windows,
    dot_aligned,
    filenameize,
    format_int_input,
    format_mass,
    loadfont,
    resolvepath,
    round_sig,
    to_si,
    unloadfont,
    validate_ce,
    validate_nn,
    validate_pi,
)
from .tip import create_tool_tip

DESCRIPTION = "Description"

root_logger = logging.getLogger()
root_logger.setLevel(logging.INFO)
handler = logging.StreamHandler(sys.stderr)
formatter = logging.Formatter("%(asctime)s - %(name)s - %(message)s")
handler.setFormatter(formatter)
root_logger.addHandler(handler)

USE_LF = "USE_LF"
USE_CV = "USE_CV"
USE_LD = "USE_LD"

FONTNAME = "Sarasa Fixed SC"
FONTSIZE = 8

CONTEXT = {
    "font.size": FONTSIZE,
    "axes.titlesize": FONTSIZE,
    "axes.labelsize": FONTSIZE,
    "axes.titlelocation": "right",
    "xtick.labelsize": FONTSIZE,
    "ytick.labelsize": FONTSIZE,
    "legend.fontsize": FONTSIZE,
    "figure.titlesize": FONTSIZE + 2,
    "lines.markersize": FONTSIZE / 4,
    "lines.linewidth": 1,
    "font.family": FONTNAME,
    "axes.labelweight": "bold",
    "xaxis.labellocation": "right",
    "yaxis.labellocation": "top",
}


THEMES = {"awlight": 1, "awdark": 0}


def grid_configure_recursive(widget, **kwargs):
    stack = list(widget.winfo_children())
    while stack:
        descendent = stack.pop()
        stack.extend(descendent.winfo_children())
        descendent.grid_configure(**kwargs)


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


class PIBS(Tk):
    def __init__(self, *args, loc, **kwargs):
        super().__init__(*args, **kwargs)
        # Only works on windows:(
        if platform.system() == "Windows":
            super().iconbitmap(default=resolvepath("ui/logo.ico"))
        super().option_add("*tearOff", False)

        font = Font(family=FONTNAME, size=FONTSIZE)
        self.option_add("*Font", font)

        dpi = self.winfo_fpixels("1i")
        # Tk was originally developed for a dpi of 72
        self.tk.call("tk", "scaling", 1.0 * dpi / 72.0)
        self.tk.call("lappend", "auto_path", resolvepath("ui/tksvg0.14"))
        self.tk.call("lappend", "auto_path", resolvepath("ui/awthemes-10.4.0"))

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
        )
        self.frame.grid(row=0, column=0, sticky="nsew")

        self.minsize(self.winfo_width(), self.winfo_height())  # set minimum size
        # self.is_fullscreen = False

        # self.bind("<F4>", lambda *_: self.toggle_fullscreen())

    # def toggle_fullscreen(self):
    #     self.wm_attributes("-fullscreen", self.is_fullscreen)
    #     self.is_fullscreen = not self.is_fullscreen

    def quit(self):
        # explicitly unload the font at the end of program.
        unloadfont(resolvepath("ui/sarasa-fixed-sc-regular.ttf"), True, True)
        super().quit()


class InteriorBallisticsFrame(LocalizedFrame):

    def __init__(self, root, menubar, dpi, default_lang, localization_dict, font, os_dark: bool):
        super().__init__(
            root, font=font, menubar=menubar, default_lang=default_lang, localization_dict=localization_dict
        )

        self.cvldlf_traces = []
        self.font = font
        self.job_queue, self.guide_job_queue, self.log_queue = Queue(), Queue(), Queue()
        self.process, self.guide_process, self.kwargs, self.gun_result = None, None, None, None
        self.dpi = dpi
        self.root = root
        self.force_update_on_theme_widget = []

        self.menubar = menubar

        data_menu = Menu(menubar)
        menubar.add_cascade(label=self.get_loc_str("dataLabel"), menu=data_menu)

        design_menu = Menu(menubar)
        menubar.add_cascade(label=self.get_loc_str("designLabel"), menu=design_menu)

        theme_menu = Menu(menubar)
        menubar.add_cascade(label=self.get_loc_str("themeLabel"), menu=theme_menu)

        debug_menu = Menu(menubar)
        menubar.add_cascade(label=self.get_loc_str("debugLabel"), menu=debug_menu)

        self.design_menu = design_menu
        # self.file_menu = file_menu
        self.data_menu = data_menu
        self.theme_menu = theme_menu
        self.debug_menu = debug_menu

        self.theme_name_var = StringVar(value="awdark" if os_dark else "awlight")
        self.debug = IntVar(value=0)

        design_menu.add_command(label=self.get_loc_str("saveLabel"), command=self.save, accelerator="Ctrl+S")
        self.root.bind("<Control-s>", lambda *_: self.save())
        self.root.bind("<Control-S>", lambda *_: self.save())

        design_menu.add_command(label=self.get_loc_str("loadLabel"), command=self.load_gun, accelerator="Ctrl+L")
        self.root.bind("<Control-l>", lambda *_: self.load_gun())
        self.root.bind("<Control-L>", lambda *_: self.load_gun())

        design_menu.add_command(
            label=self.get_loc_str("loadPresetLabel"),
            command=lambda: self.load_gun(initial_dir=resolvepath("examples")),
        )

        design_menu.add_command(label=self.get_loc_str("resetLabel"), command=self.on_reset, accelerator="Ctrl+N")
        self.root.bind("<Control-N>", lambda *_: self.on_reset())
        self.root.bind("<Control-n>", lambda *_: self.on_reset())

        design_menu.add_command(label=self.get_loc_str("calcLabel"), command=self.on_calculate, accelerator="Ctrl+R")
        self.root.bind("<Control-R>", lambda *_: self.on_calculate())
        self.root.bind("<Control-r>", lambda *_: self.on_calculate())

        design_menu.add_command(label=self.get_loc_str("guideLabel"), command=self.on_guide, accelerator="Ctrl+G")
        self.root.bind("<Control-G>", lambda *_: self.on_guide())
        self.root.bind("<Control-g>", lambda *_: self.on_guide())

        data_menu.add_command(label=self.get_loc_str("exportMain"), command=lambda: self.export_graph(save="main"))
        data_menu.add_command(label=self.get_loc_str("exportAux"), command=lambda: self.export_graph(save="aux"))
        data_menu.add_command(label=self.get_loc_str("exportGeom"), command=lambda: self.export_graph(save="geom"))
        data_menu.add_command(label=self.get_loc_str("exportGuide"), command=lambda: self.export_graph(save="guide"))
        data_menu.add_command(label=self.get_loc_str("exportLabel"), command=self.export_table)

        data_menu.add_command(label=self.get_loc_str("reloadPropellant"), command=self.load_propellant)

        for themeName in THEMES.keys():
            theme_menu.add_radiobutton(
                label=themeName, variable=self.theme_name_var, value=themeName, command=self.use_theme
            )

        debug_menu.add_checkbutton(label=self.get_loc_str("enableLabel"), variable=self.debug, onvalue=1, offvalue=0)

        self.prop, self.gun, self.guide = None, None, None

        self.columnconfigure(1, weight=1)
        self.rowconfigure(0, weight=1)
        self.rowconfigure(1, weight=2)

        ## nameplate
        name_frm = self.add_localized_label_frame(self, label_loc_key="nameFrm")
        name_frm.grid(row=0, column=0, sticky="nsew", columnspan=2)
        name_frm.columnconfigure(0, weight=1)
        name_frm.rowconfigure(1, weight=1)

        self.name = StringVar(self)
        name_plate = ttk.Entry(
            name_frm, textvariable=self.name, state="disabled", justify="left", font=(FONTNAME, FONTSIZE + 2)
        )
        name_plate.grid(row=0, column=0, sticky="nsew", padx=2, pady=2, columnspan=2)

        desc_scroll = ttk.Scrollbar(name_frm, orient="vertical")
        desc_scroll.grid(row=1, column=1, sticky="nsew")
        self.description = Text(
            name_frm, wrap="word", height=5, width=80, yscrollcommand=desc_scroll.set, font=(FONTNAME, FONTSIZE)
        )
        self.description.grid(row=1, column=0, sticky="nsew")

        desc_scroll.config(command=self.description.yview)
        self.force_update_on_theme_widget.append(self.description)

        ## left frame
        info_frame = Frame(self)
        info_frame.grid(row=1, column=0, sticky="nsew")
        info_frame.columnconfigure(0, weight=1)
        info_frame.rowconfigure(0, weight=1)

        par_frm = self.add_localized_label_frame(info_frame, label_loc_key="parFrmLabel")
        par_frm.grid(row=0, column=0, sticky="nsew")
        par_frm.columnconfigure(0, weight=1)

        i = 0
        self.ammo, i = (
            self.add_localized_12_display(parent=par_frm, row=i, label_loc_key="ammoLabel", unit_text=""),
            i + 2,
        )
        self.pp, i = (
            self.add_localized_122_display(parent=par_frm, row=i, label_loc_key="ppLabel", tooltip_loc_key="ppText"),
            i + 3,
        )
        self.bop, i = self.add_localized_12_display(parent=par_frm, row=i, label_loc_key="bopLabel"), i + 2
        self.lx, i = (
            self.add_localized_122_display(parent=par_frm, row=i, label_loc_key="lxLabel", tooltip_loc_key="calLxText"),
            i + 3,
        )
        self.mv, i = self.add_localized_12_display(parent=par_frm, row=i, label_loc_key="mvLabel"), i + 2
        self.va, i = (
            self.add_localized_12_display(parent=par_frm, row=i, label_loc_key="vaLabel", tooltip_loc_key="vinfText"),
            i + 2,
        )
        self.te, i = (
            self.add_localized_12_display(parent=par_frm, row=i, label_loc_key="teffLabel", tooltip_loc_key="teffText"),
            i + 2,
        )
        self.be, i = (
            self.add_localized_12_display(parent=par_frm, row=i, label_loc_key="beffLabel", tooltip_loc_key="beffText"),
            i + 2,
        )
        self.pe, i = (
            self.add_localized_12_display(parent=par_frm, row=i, label_loc_key="peffLabel", tooltip_loc_key="peffText"),
            i + 2,
        )
        self.pa, i = self.add_localized_12_display(parent=par_frm, row=i, label_loc_key="paLabel"), i + 2
        self.gm, i = self.add_localized_12_display(parent=par_frm, row=i, label_loc_key="gmLabel"), i + 2

        self.sj, i = self.add_localized_12_display(parent=par_frm, row=i, label_loc_key="sjLabel"), i + 2

        ## center frame
        self.tab_parent = ttk.Notebook(self, padding=0)
        self.tab_parent.grid(row=1, column=1, sticky="nsew")

        self.tab_parent.columnconfigure(0, weight=1)
        self.tab_parent.rowconfigure(0, weight=1)

        self.plot_tab = Frame(self.tab_parent)
        self.plot_tab.grid(row=0, column=0, stick="nsew")
        self.plot_tab.rowconfigure(0, weight=2)
        self.plot_tab.rowconfigure(1, weight=1)
        self.plot_tab.columnconfigure(0, weight=1)

        self.table_tab = Frame(self.tab_parent)
        self.table_tab.grid(row=0, column=0, sticky="nsew")
        self.table_tab.rowconfigure(0, weight=1)
        self.table_tab.columnconfigure(0, weight=1)

        self.guide_tab = Frame(self.tab_parent)
        self.guide_tab.grid(row=0, column=0, sticky="nsew")
        self.guide_tab.rowconfigure(0, weight=1)
        self.guide_tab.columnconfigure(0, weight=1)

        self.errorTab = Frame(self.tab_parent)
        self.errorTab.grid(row=0, column=0, sticky="nsew")
        self.errorTab.rowconfigure(0, weight=1)
        self.errorTab.columnconfigure(0, weight=1)

        self.tab_parent.add(self.plot_tab, text=self.get_loc_str("plotTab"))
        self.tab_parent.add(self.table_tab, text=self.get_loc_str("tableTab"))
        self.tab_parent.add(self.guide_tab, text=self.get_loc_str("guideTab"))
        self.tab_parent.add(self.errorTab, text=self.get_loc_str("errorTab"))

        self.tab_parent.enable_traversal()

        ## table frame
        tbl_frm = self.add_localized_label_frame(self.table_tab, label_loc_key="tblFrmLabel")
        tbl_frm.grid(row=0, column=0, sticky="nsew")

        tbl_frm.columnconfigure(0, weight=1)
        tbl_frm.rowconfigure(0, weight=1)

        tbl_place_frm = Frame(tbl_frm)
        tbl_place_frm.grid(row=0, column=0, sticky="nsew")

        vertscroll = ttk.Scrollbar(tbl_frm, orient="vertical")  # create a scrollbar
        vertscroll.grid(row=0, rowspan=2, column=1, sticky="nsew")

        horzscroll = ttk.Scrollbar(tbl_frm, orient="horizontal")
        horzscroll.grid(row=1, column=0, sticky="nsew")

        self.tv = ttk.Treeview(
            tbl_place_frm, selectmode="browse", yscrollcommand=vertscroll.set, xscrollcommand=horzscroll.set
        )  # this set the nbr. of values
        self.tv.place(relwidth=1, relheight=1)

        vertscroll.configure(command=self.tv.yview)  # make it vertical
        horzscroll.configure(command=self.tv.xview)

        ## error frame
        error_frm = self.add_localized_label_frame(self.errorTab, label_loc_key="errFrmLabel")
        error_frm.grid(row=0, column=0, columnspan=3, sticky="nsew")
        error_frm.columnconfigure(0, weight=1)
        error_frm.rowconfigure(0, weight=1)

        err_scroll = ttk.Scrollbar(error_frm, orient="vertical")
        err_scroll.grid(row=0, column=1, sticky="nsew")
        self.error_text = Text(
            error_frm, yscrollcommand=err_scroll.set, wrap="word", height=6, width=0, font=(FONTNAME, FONTSIZE)
        )
        self.error_text.grid(row=0, column=0, sticky="nsew")

        self.error_text.tag_configure(str(logging.DEBUG), foreground="tan")
        self.error_text.tag_configure(str(logging.WARNING), foreground="orange")
        self.error_text.tag_configure(str(logging.ERROR), foreground="orangered")
        self.error_text.tag_configure(str(logging.CRITICAL), foreground="red")

        err_scroll.config(command=self.error_text.yview)
        self.force_update_on_theme_widget.append(self.error_text)

        plot_frm = self.add_localized_label_frame(
            self.plot_tab, label_loc_key="plotFrmLabel", tooltip_loc_key="plotText"
        )
        plot_frm.grid(row=0, column=0, sticky="nsew")

        for i in range(5):
            plot_frm.columnconfigure(i, weight=1)

        j = 1
        self.plot_avg_p, self.plot_base_p, self.plot_breech_p, self.plot_stag_p = (
            self.add_localized_label_check(parent=plot_frm, label_loc_key=label, desc_label_key=None, row=j, col=k)
            for k, label in enumerate(("plotAvgP", "plotBaseP", "plotBreechP", "plotStagP"))
        )

        j += 1
        self.plot_vel, self.plot_nozzle_v, self.plot_burnup, self.plot_eta = (
            self.add_localized_label_check(parent=plot_frm, label_loc_key=label, desc_label_key=None, row=j, col=k)
            for k, label in enumerate(("plotVel", "plotNozzleV", "plotBurnup", "plotEta"))
        )

        plot_frm.columnconfigure(0, weight=1)
        plot_frm.rowconfigure(0, weight=1)

        plot_place_frm = Frame(plot_frm)
        plot_place_frm.grid(row=0, column=0, padx=2, pady=2, sticky="nsew", columnspan=5)

        with mpl.rc_context(CONTEXT):
            fig = Figure(dpi=96, layout="constrained")
            axes = fig.add_subplot(111)

            ax = axes
            ax_p, ax_v = ax.twinx(), ax.twinx()

            ax.yaxis.tick_right()
            ax_v.yaxis.tick_left()

            ax.set_xlabel(" ")
            # noinspection PyTypeChecker
            ax_p.spines.right.set_position(("data", 0.5))
            ax_p.yaxis.set_ticks(ax_p.get_yticks()[1:-1:])

            self.ax, self.ax_p, self.ax_v, self.fig = ax, ax_p, ax_v, fig

            self.plt_canvas = FigureCanvasTkAgg(fig, master=plot_place_frm)
            self.plt_canvas.draw_idle()
            self.plt_canvas.get_tk_widget().place(relheight=1, relwidth=1)

        aux_frm = self.add_localized_label_frame(self.plot_tab, label_loc_key="auxFrmLabel", tooltip_loc_key="auxText")
        aux_frm.grid(row=1, column=0, sticky="nsew")

        for i in range(2):
            aux_frm.columnconfigure(i, weight=1)

        j = 1
        k = 0
        self.trace_hull, k = (
            self.add_localized_label_check(
                parent=aux_frm, row=j, col=k, label_loc_key="traceHull", default=True, desc_label_key=None
            ),
            k + 1,
        )

        self.trace_press, k = (
            self.add_localized_label_check(
                parent=aux_frm, row=j, col=k, label_loc_key="tracePress", desc_label_key=None
            ),
            k + 1,
        )

        aux_frm.columnconfigure(0, weight=1)
        aux_frm.rowconfigure(0, weight=1)

        aux_place_frm = Frame(aux_frm)
        aux_place_frm.grid(row=0, column=0, padx=2, pady=2, sticky="nsew", columnspan=2)

        with mpl.rc_context(CONTEXT):
            fig = Figure(dpi=96, layout="constrained")

            self.aux_ax = fig.add_subplot(111)
            self.aux_ax_h = self.aux_ax.twinx()
            self.aux_fig = fig
            self.aux_ax_h.yaxis.tick_right()
            self.aux_canvas = FigureCanvasTkAgg(fig, master=aux_place_frm)
            self.aux_canvas.draw_idle()
            self.aux_canvas.get_tk_widget().place(relwidth=1, relheight=1)

        ## right frame
        right_frm = Frame(self)
        right_frm.grid(row=0, column=4, rowspan=2, sticky="nsew")
        right_frm.columnconfigure(0, weight=1)
        right_frm.rowconfigure(0, weight=1)

        validation_nn = self.register(validate_nn)
        validation_pi = self.register(validate_pi)
        validation_ce = self.register(validate_ce)

        i = 1
        env_frm = self.add_localized_label_frame(right_frm, label_loc_key="envFrmLabel")
        env_frm.grid(row=i, column=0, sticky="nsew")
        i += 1
        env_frm.columnconfigure(0, weight=1)

        j = 0
        self.in_atmos, j = (
            self.add_localized_label_check(
                parent=env_frm, label_loc_key="atmosLabel", desc_label_key="atmosLabel", row=j, columnspan=3
            ),
            j + 1,
        )

        self.amb_p, j = (
            self.add_localized_3_input(
                parent=env_frm,
                row=j,
                label_loc_key="ambPresLabel",
                unit_text="kPa",
                default="101.325",
                validation=validation_nn,
                dtype=float,
            ),
            j + 1,
        )
        self.amb_rho, j = (
            self.add_localized_3_input(
                parent=env_frm,
                row=j,
                label_loc_key="ambRhoLabel",
                unit_text="kg/m³",
                default="1.204",
                validation=validation_nn,
                dtype=float,
            ),
            j + 1,
        )

        self.amb_gamma, j = (
            self.add_localized_3_input(
                parent=env_frm,
                row=j,
                label_loc_key="ambGamLabel",
                default="1.400",
                validation=validation_nn,
                dtype=float,
            ),
            j + 1,
        )

        sol_frm = self.add_localized_label_frame(right_frm, label_loc_key="solFrmLabel")
        sol_frm.grid(row=i, column=0, sticky="nsew")
        i += 1
        sol_frm.columnconfigure(0, weight=1)

        self.drop_soln = self.add_localized_dropdown(
            parent=sol_frm,
            str_obj_dict={soln: soln for soln in (SOL_LAGRANGE, SOL_PIDDUCK, SOL_MAMONTOV)},
            desc_label_key="solFrmLabel",
        )
        self.drop_soln.grid(row=0, column=0, columnspan=2, sticky="nsew", padx=2, pady=2)

        op_frm = self.add_localized_label_frame(right_frm, label_loc_key="opFrmLabel")
        op_frm.grid(row=4, column=0, sticky="nsew")
        op_frm.columnconfigure(0, weight=1)

        self.solve_W_Lg, i = (
            self.add_localized_label_check(
                parent=op_frm,
                row=i,
                columnspan=3,
                default=False,
                label_loc_key="consButton",
                desc_label_key="consButton",
                tooltip_loc_key="useConsText",
            ),
            i + 1,
        )

        cons_frm = self.add_localized_label_frame(
            op_frm, label_loc_key="consFrmLabel"
        )  # style="SubLabelFrame.TLabelframe"
        cons_frm.grid(row=i, column=0, columnspan=3, sticky="nsew", padx=2, pady=2)
        i += 1
        cons_frm.columnconfigure(0, weight=1)

        j = 0
        self.lock_Lg, j = (
            self.add_localized_label_check(
                parent=cons_frm,
                row=j,
                columnspan=3,
                default=False,
                label_loc_key="lockButton",
                desc_label_key="lockButton",
                tooltip_loc_key="lockText",
            ),
            j + 1,
        )

        self.opt_lf, j = (
            self.add_localized_label_check(
                parent=cons_frm,
                row=j,
                columnspan=3,
                default=False,
                desc_label_key="minTVButton",
                label_loc_key="minTVButton",
                tooltip_loc_key="optLFText",
            ),
            j + 1,
        )

        self.v_tgt, j = (
            self.add_localized_3_input(
                parent=cons_frm,
                row=j,
                label_loc_key="vTgtLabel",
                unit_text="m/s",
                default="1000.0",
                validation=validation_nn,
                dtype=float,
            ),
            j + 1,
        )

        self.p_tgt, j = (
            self.add_localized_3_input(
                parent=cons_frm,
                row=j,
                label_loc_key="pTgtLabel",
                unit_text="MPa",
                default="350.0",
                validation=validation_nn,
                tooltip_loc_key="pTgtText",
                dtype=float,
            ),
            j + 1,
        )

        self.p_control = self.add_localized_dropdown(parent=cons_frm, desc_label_key="Pressure Constraint")
        self.p_control.grid(row=j, column=0, columnspan=3, sticky="nsew", padx=2, pady=2)
        j += 1

        self.min_web, j = (
            self.add_localized_3_input(
                parent=cons_frm,
                row=j,
                label_loc_key="minWebLabel",
                unit_text="μm",
                default="1.0",
                validation=validation_nn,
                color="red",
                dtype=float,
            ),
            j + 1,
        )
        self.lg_max, j = (
            self.add_localized_3_input(
                parent=cons_frm,
                row=j,
                label_loc_key="maxLgLabel",
                unit_text="m",
                default="10.0",
                validation=validation_nn,
                color="red",
                dtype=float,
            ),
            j + 1,
        )

        sample_frm = self.add_localized_label_frame(
            op_frm, label_loc_key="sampleFrmLabel", style="SubLabelFrame.TLabelframe", tooltip_loc_key="sampText"
        )
        sample_frm.grid(row=i, column=0, columnspan=2, sticky="nsew", padx=2, pady=2)
        sample_frm.columnconfigure(0, weight=1)
        i += 1

        j = 0

        self.drop_domain = self.add_localized_dropdown(
            parent=sample_frm,
            str_obj_dict={domain: domain for domain in (DOMAIN_TIME, DOMAIN_LEN)},
            desc_label_key="sampleFrmLabel",
        )
        self.drop_domain.grid(row=j, column=0, columnspan=2, sticky="nsew", padx=2, pady=2)
        j += 1

        self.step, j = (
            self.add_localized_2_input(
                parent=sample_frm,
                row=j,
                label_loc_key="stepLabel",
                default="33",
                validation=validation_nn,
                formatter=format_int_input,
            ),
            j + 1,
        )

        self.acc_exp, i = (
            self.add_localized_2_input(
                parent=op_frm,
                row=i,
                label_loc_key="-log10(ε)",
                default="3",
                validation=validation_pi,
                formatter=format_int_input,
                color="red",
                tooltip_loc_key="tolText",
            ),
            i + 1,
        )

        self.max_iter, i = (
            self.add_localized_2_input(
                parent=op_frm,
                row=i,
                label_loc_key="maxIterLabel",
                default="10",
                validation=validation_pi,
                formatter=format_int_input,
                color="red",
            ),
            i + 1,
        )

        self.max_guess, i = (
            self.add_localized_2_input(
                parent=op_frm,
                row=i,
                label_loc_key="maxGuessLabel",
                default="25",
                validation=validation_pi,
                formatter=format_int_input,
                color="red",
            ),
            i + 1,
        )

        self.calc_button = ttk.Button(op_frm, text=self.get_loc_str("calcLabel"), command=self.on_calculate)
        self.calc_button.grid(row=i, column=0, columnspan=3, sticky="nsew", padx=2, pady=2)

        op_frm.rowconfigure(i, weight=1)
        i += 1

        self.calc_button_tip = StringVar(value=self.get_loc_str("calcButtonText"))
        create_tool_tip(self.calc_button, self.calc_button_tip, font=self.font)

        ctrl_frm = self.add_localized_label_frame(right_frm, label_loc_key="guideCtrlFrmLabel")
        ctrl_frm.grid(row=i, column=0, columnspan=3, sticky="nsew")
        i += 1
        ctrl_frm.columnconfigure(0, weight=1)

        (
            self.guide_min_lf,
            self.guide_max_lf,
            self.guide_step_lf,
            self.guide_min_cmr,
            self.guide_max_cmr,
            self.guide_step_cmr,
        ) = (
            self.add_localized_3_input(
                ctrl_frm,
                label_loc_key=locKey,
                default=default,
                unit_text=unit,
                validation=validation,
                row=j,
                dtype=float,
            )
            for j, (locKey, default, unit, validation) in enumerate(
                (
                    ("minLFLabel", "10.0", "%", validation_ce),
                    ("maxLFLabel", "90.0", "%", validation_ce),
                    ("stepLFLabel", "5.0", "%", validation_ce),
                    ("minCMRLabel", "0.25", "", validation_nn),
                    ("maxCMRLabel", "1.75", "", validation_nn),
                    ("stepCMRLabel", "0.05", "", validation_nn),
                )
            )
        )
        self.guide_button = ttk.Button(ctrl_frm, text=self.get_loc_str("guideLabel"), command=self.on_guide)

        self.guide_button.grid(row=6, column=0, columnspan=3, sticky="nsew", padx=2, pady=2)

        left_frame = Frame(self)
        left_frame.grid(row=0, column=2, rowspan=2, sticky="nsew")
        left_frame.columnconfigure(0, weight=1)
        left_frame.rowconfigure(0, weight=1)

        ## spec frame
        specs_frame = self.add_localized_label_frame(left_frame, label_loc_key="specFrmLabel")
        specs_frame.grid(row=0, column=0, sticky="nsew")
        specs_frame.columnconfigure(0, weight=1)

        validation_nn, validation_ce = self.register(validate_nn), self.register(validate_ce)

        i = 0
        self.type_optn = self.add_localized_dropdown(
            parent=specs_frame,
            str_obj_dict={gun_type: gun_type for gun_type in (CONVENTIONAL, RECOILLESS)},
            desc_label_key="typeLabel",
        )
        self.type_optn.grid(row=i, column=0, sticky="nsew", padx=2, pady=2, columnspan=3)
        i += 1

        self.cal_mm, i = (
            self.add_localized_3_input(
                parent=specs_frame,
                row=i,
                label_loc_key="calLabel",
                unit_text="mm",
                default="50.0",
                validation=validation_nn,
                dtype=float,
            ),
            i + 1,
        )

        self.tbl_mm, i = (
            self.add_localized_3_input(
                parent=specs_frame,
                row=i,
                label_loc_key="tblLabel",
                unit_text="mm",
                default="3500.0",
                validation=validation_nn,
                dtype=float,
            ),
            i + 1,
        )

        self.sht_kg, i = (
            self.add_localized_3_input(
                parent=specs_frame,
                row=i,
                label_loc_key="shtLabel",
                unit_text="kg",
                default="2.0",
                validation=validation_nn,
                dtype=float,
            ),
            i + 1,
        )

        self.chg_kg, i = (
            self.add_localized_3_input(
                parent=specs_frame,
                row=i,
                label_loc_key="chgLabel",
                unit_text="kg",
                default="0.5",
                validation=validation_nn,
                tooltip_loc_key="chgText",
                dtype=float,
            ),
            i + 1,
        )

        self.use_cv = self.add_localized_dropdown(
            parent=specs_frame,
            str_obj_dict={USE_CV: USE_CV, USE_LD: USE_LD, USE_LF: USE_LF},
            desc_label_key="cvlfLabel",
        )
        self.use_cv.grid(row=i, column=0, columnspan=3, sticky="nsew", padx=2, pady=2)
        i += 1

        self.cv_L, i = (
            self.add_localized_3_input(
                parent=specs_frame,
                row=i,
                label_loc_key="cvLabel",
                unit_text="L",
                default="1.0",
                validation=validation_nn,
                dtype=float,
            ),
            i + 1,
        )

        self.ld, i = (
            self.add_localized_3_input(
                parent=specs_frame,
                row=i,
                label_loc_key="ldLabel",
                unit_text="kg/m³",
                default="0.0",
                validation=validation_nn,
                dtype=float,
            ),
            i + 1,
        )

        self.lf, i = (
            self.add_localized_3_input(
                parent=specs_frame,
                row=i,
                label_loc_key="ldfLabel",
                unit_text="%",
                default="0.0",
                validation=validation_ce,
                tooltip_loc_key="ldfText",
                dtype=float,
            ),
            i + 1,
        )

        self.clr, i = (
            self.add_localized_3_input(
                parent=specs_frame,
                row=i,
                label_loc_key="clrLabel",
                unit_text="x",
                default="1.5",
                validation=validation_nn,
                tooltip_loc_key="clrText",
                dtype=float,
            ),
            i + 1,
        )

        self.dgc, i = (
            self.add_localized_3_input(
                parent=specs_frame,
                row=i,
                label_loc_key="dgcLabel",
                unit_text="%",
                default="3.0",
                validation=validation_ce,
                tooltip_loc_key="dgcText",
                dtype=float,
            ),
            i + 1,
        )

        self.stp_MPa, i = (
            self.add_localized_3_input(
                parent=specs_frame,
                row=i,
                label_loc_key="stpLabel",
                unit_text="MPa",
                default="30.0",
                validation=validation_nn,
                tooltip_loc_key="stpText",
                dtype=float,
            ),
            i + 1,
        )

        self.nozz_exp, i = (
            self.add_localized_3_input(
                parent=specs_frame,
                row=i,
                label_loc_key="nozzExpLabel",
                unit_text="x",
                default="4.0",
                validation=validation_nn,
                tooltip_loc_key="nozzExpText",
                dtype=float,
            ),
            i + 1,
        )

        self.nozz_eff, i = (
            self.add_localized_3_input(
                parent=specs_frame,
                row=i,
                label_loc_key="nozzEffLabel",
                unit_text="%",
                default="92.0",
                validation=validation_ce,
                tooltip_loc_key="nozzEffText",
                dtype=float,
            ),
            i + 1,
        )

        self.use_material, i = (
            self.add_localized_label_check(
                specs_frame,
                row=i,
                label_loc_key="useMaterialLabel",
                desc_label_key="useMaterialLabel",
                default=False,
                columnspan=3,
            ),
            i + 1,
        )

        material_frame = self.add_localized_label_frame(specs_frame, label_loc_key="matFrmLabel")
        material_frame.grid(row=i, column=0, sticky="nsew", columnspan=3, padx=2, pady=2)
        material_frame.columnconfigure(0, weight=1)
        i += 1

        j = 0

        self.material_density, j = (
            self.add_localized_3_input(
                material_frame,
                row=j,
                label_loc_key="matDensityLabel",
                unit_text="kg/m³",
                default="7850.0",
                validation=validation_nn,
                dtype=float,
            ),
            j + 1,
        )
        self.material_yield, j = (
            self.add_localized_3_input(
                material_frame,
                row=j,
                label_loc_key="matYieldLabel",
                unit_text="MPa",
                default="1000.0",
                validation=validation_nn,
                dtype=float,
            ),
            j + 1,
        )

        self.material_ssf, j = (
            self.add_localized_2_input(
                parent=material_frame,
                row=j,
                label_loc_key="sffLabel",
                default="1.35",
                validation=validation_nn,
                dtype=float,
            ),
            j + 1,
        )

        self.material_is_af, i = (
            self.add_localized_label_check(
                parent=material_frame, label_loc_key="afLabel", desc_label_key="afLabel", row=j, columnspan=2
            ),
            i + 1,
        )

        ## grain frame

        mid_frame = Frame(self)

        mid_frame.grid(row=0, column=3, rowspan=2, sticky="nsew")
        mid_frame.columnconfigure(0, weight=1)

        i = 0
        mid_frame.rowconfigure(i, weight=1)
        ### propellant frame
        prop_frame = self.add_localized_label_frame(
            mid_frame, label_loc_key="propFrmLabel", tooltip_loc_key="specsText"
        )
        prop_frame.grid(row=i, column=0, columnspan=3, sticky="nsew")
        i += 1
        prop_frame.rowconfigure(1, weight=1)
        prop_frame.columnconfigure(0, weight=1)
        j = 0
        self.drop_prop = self.add_localized_dropdown(
            parent=prop_frame,
            str_obj_dict=Composition.read_file(resolvepath("ballistics/resource/propellants.csv")),
            desc_label_key="propFrmLabel",
        )
        self.drop_prop.grid(row=j, column=0, columnspan=2, sticky="nsew")
        j += 1

        spec_scroll = ttk.Scrollbar(prop_frame, orient="vertical")
        spec_h_scroll = ttk.Scrollbar(prop_frame, orient="horizontal")
        self.specs = Text(
            prop_frame,
            wrap="word",
            height=0,
            width=0,
            yscrollcommand=spec_scroll.set,
            xscrollcommand=spec_h_scroll.set,
            font=(FONTNAME, FONTSIZE),
        )
        self.force_update_on_theme_widget.append(self.specs)
        spec_scroll.config(command=self.specs.yview)
        spec_h_scroll.config(command=self.specs.xview)

        self.specs.grid(row=j, column=0, sticky="nsew")
        spec_scroll.grid(row=j, rowspan=2, column=1, sticky="nsew")
        j += 1
        spec_h_scroll.grid(row=j, column=0, sticky="nsew")
        j += 1

        self.use_combustible, j = (
            self.add_localized_label_check(
                prop_frame,
                label_loc_key="combustibleLabel",
                tooltip_loc_key="combustibleText",
                default=False,
                row=j,
                columnspan=2,
            ),
            j + 1,
        )

        combustible_frame = Frame(prop_frame)
        combustible_frame.grid(row=j, column=0, columnspan=2, sticky="nsew")
        j += 1

        k = 0
        self.combustible_mass_kg, k = (
            self.add_localized_3_input(
                combustible_frame,
                label_loc_key="combustibleMassLabel",
                default="0.0",
                unit_text="kg",
                desc_label_key="ωʹ",
                dtype=float,
                validation=validation_nn,
                row=k,
            ),
            k + 1,
        )

        self.combustible_force_kJ__kg, k = (
            self.add_localized_3_input(
                combustible_frame,
                label_loc_key="combustibleForceLabel",
                default="750.0",
                unit_text="kJ/kg",
                desc_label_key="fʹ",
                dtype=float,
                validation=validation_nn,
                row=k,
            ),
            k + 1,
        )
        force_fudge_frame = Frame(prop_frame)
        force_fudge_frame.grid(row=j, column=0, columnspan=2, sticky="nsew")
        j += 1

        force_fudge_frame.columnconfigure(0, weight=1)
        force_fudge_frame.rowconfigure(0, weight=1)

        self.force_fudge = self.add_localized_3_input(
            force_fudge_frame,
            label_loc_key="forceFudgeLabel",
            tooltip_loc_key="forceFudgeText",
            default="100.0",
            unit_text="%",
            dtype=float,
            validation=validation_nn,
        )
        ###

        mid_frame.rowconfigure(i, weight=2)

        self.grain_frm = self.add_localized_label_frame(mid_frame, label_loc_key="grainFrmLabel")
        self.grain_frm.grid(row=i, column=0, sticky="nsew")
        i += 1
        self.grain_frm.columnconfigure(0, weight=1)
        self.grain_frm.rowconfigure(0, weight=1)

        j = 0
        self.geom_plot_frm = self.add_localized_label_frame(
            self.grain_frm, label_loc_key="σ(Z)", style="SubLabelFrame.TLabelframe", tooltip_loc_key="geomPlotText"
        )
        self.geom_plot_frm.grid(row=j, column=0, columnspan=3, sticky="nsew", padx=2, pady=2)
        j += 1

        ## geom plot
        with mpl.rc_context(CONTEXT):
            fig = Figure(dpi=96, layout="constrained")
            self.geom_fig = fig
            self.geom_ax = fig.add_subplot(111)

            self.geom_canvas = FigureCanvasTkAgg(fig, master=self.geom_plot_frm)
            self.geom_canvas.draw_idle()
            self.geom_canvas.get_tk_widget().place(relheight=1, relwidth=1)

        self.main_geom = self.add_localized_dropdown(
            parent=self.grain_frm, str_obj_dict=Geometry.get_desc_geometry_dict(), desc_label_key="Grain Geometry"
        )
        self.main_geom.grid(row=j, column=0, columnspan=3, sticky="nsew", padx=2, pady=2)
        j += 1

        self.web_mm, j = (
            self.add_localized_3_input(
                parent=self.grain_frm,
                row=j,
                desc_label_key="Web",
                unit_text="mm",
                default="1.0",
                validation=validation_nn,
                dtype=float,
            ),
            j + 1,
        )

        self.grain_r1, j = (
            self.add_localized_3_input(
                parent=self.grain_frm,
                row=j,
                desc_label_key="1/α",
                unit_text="x",
                default="1.0",
                validation=validation_nn,
                tooltip_loc_key="",
                dtype=float,
            ),
            j + 1,
        )

        self.grain_r2, j = (
            self.add_localized_3_input(
                parent=self.grain_frm,
                row=j,
                desc_label_key="1/β",
                unit_text="x",
                default="10.0",
                validation=validation_nn,
                dtype=float,
            ),
            j + 1,
        )

        self.use_aux_grain, j = (
            self.add_localized_label_check(
                parent=self.grain_frm,
                row=j,
                label_loc_key="useAuxGrainLabel",
                desc_label_key="useAuxGrainLabel",
                columnspan=3,
                default=False,
            ),
            j + 1,
        )

        self.aux_grain_frm = self.add_localized_label_frame(self.grain_frm, label_loc_key="auxGrainFrmLabel")
        self.aux_grain_frm.grid(row=j, column=0, columnspan=3, sticky="nsew", padx=2, pady=2)
        j += 1
        self.aux_grain_frm.columnconfigure(1, weight=1)

        k = 0
        self.aux_mass_ratio, k = (
            self.add_localized_3_input(
                parent=self.aux_grain_frm,
                row=k,
                label_loc_key="auxMassRatio",
                default="1.0",
                unit_text="x",
                validation=validation_nn,
                dtype=float,
            ),
            k + 1,
        )

        self.aux_geom = self.add_localized_dropdown(
            parent=self.aux_grain_frm,
            str_obj_dict=Geometry.get_desc_geometry_dict(),
            desc_label_key="Auxiliary Grain Geometry",
        )
        self.aux_geom.grid(row=k, column=0, columnspan=3, sticky="nsew", padx=2, pady=2)
        k += 1

        self.aux_web_ratio, k = (
            self.add_localized_3_input(
                parent=self.aux_grain_frm,
                row=k,
                default="1.0",
                unit_text="x",
                label_loc_key="auxWebRatio",
                validation=validation_nn,
                dtype=float,
            ),
            k + 1,
        )

        self.aux_grain_r1, k = (
            self.add_localized_3_input(
                parent=self.aux_grain_frm,
                row=k,
                unit_text="x",
                default="1.0",
                desc_label_key="Auxiliary 1/α",
                validation=validation_nn,
                tooltip_loc_key="",
                dtype=float,
            ),
            k + 1,
        )

        self.aux_grain_r2, k = (
            self.add_localized_3_input(
                parent=self.aux_grain_frm,
                row=k,
                unit_text="x",
                default="10.0",
                desc_label_key="Auxiliary 1/β",
                validation=validation_nn,
                dtype=float,
            ),
            k + 1,
        )

        self.swap_button = ttk.Button(self.grain_frm, text=self.get_loc_str("swapLabel"), command=self.swap)
        self.swap_button.grid(row=j, column=0, columnspan=3, sticky="nsew", padx=2, pady=2)
        # self.locs.append(self.swap_button)
        j += 1

        ## guide plot
        plot_frm = self.add_localized_label_frame(self.guide_tab, label_loc_key="guideFrmLabel")
        plot_frm.grid(row=0, column=0, sticky="nsew")

        with mpl.rc_context(CONTEXT):
            fig = Figure(dpi=96, layout="constrained")
            self.guide_fig = fig
            self.guide_ax = fig.add_subplot(111)
            self.guide_canvas = FigureCanvasTkAgg(fig, master=plot_frm)
            self.guide_canvas.draw_idle()
            self.guide_canvas.get_tk_widget().place(relheight=1, relwidth=1)

        control_frm = Frame(self.guide_tab)
        control_frm.grid(row=1, column=0, sticky="nsew")

        self.guide_index = IntVar(value=3)

        for i in range(3):
            control_frm.columnconfigure(i, weight=1)
        self.guide_plot_travel = ttk.Radiobutton(
            control_frm, text=self.get_loc_str("guidePlotTravel"), variable=self.guide_index, value=3
        )
        self.guide_plot_travel.grid(row=0, column=0, sticky="nsew")

        self.guide_plot_volume = ttk.Radiobutton(
            control_frm, text=self.get_loc_str("guidePlotVolume"), variable=self.guide_index, value=4
        )
        self.guide_plot_volume.grid(row=0, column=1, sticky="nsew")

        self.guide_plot_burnout = ttk.Radiobutton(
            control_frm, text=self.get_loc_str("guidePlotBurnout"), variable=self.guide_index, value=5
        )
        self.guide_plot_burnout.grid(row=0, column=2, sticky="nsew")

        self.solve_W_Lg.trace_add("write", self.ctrl_callback)
        self.lock_Lg.trace_add("write", self.ctrl_callback)
        self.opt_lf.trace_add("write", self.ctrl_callback)
        self.use_aux_grain.trace_add("write", self.ctrl_callback)

        self.guide_index.trace_add("write", callback=self.guide_callback)

        self.in_atmos.trace_add("write", self.amb_callback)

        self.drop_prop.trace_add("write", self.update_spec)

        self.drop_prop.trace_add("write", self.cvldlf_consistency_callback)

        self.main_geom.trace_add("write", self.update_geom)
        self.aux_geom.trace_add("write", self.update_geom)

        self.type_optn.trace_add("write", self.type_callback)

        self.use_material.trace_add("write", self.material_callback)

        self.use_cv.trace_add("write", self.cvldlf_callback)

        self.add_cvldlf_consistency_traces()

        for entry in (
            self.main_geom,
            self.aux_geom,
            self.use_aux_grain,
            self.drop_prop,
            self.grain_r1,
            self.grain_r2,
            self.web_mm,
            self.aux_mass_ratio,
            self.aux_web_ratio,
            self.aux_grain_r1,
            self.aux_grain_r2,
            self.use_combustible,
            self.combustible_force_kJ__kg,
            self.combustible_mass_kg,
            self.chg_kg,
            self.force_fudge,
        ):

            entry.var.trace_add("write", self.propellant_callback)

        for check in (
            self.plot_avg_p,
            self.plot_base_p,
            self.plot_breech_p,
            self.plot_stag_p,
            self.plot_vel,
            self.plot_nozzle_v,
            self.plot_burnup,
            self.plot_eta,
        ):
            check.trace_add("write", self.update_fig_plot)

        for check in (self.trace_hull, self.trace_press):
            check.trace_add("write", self.update_aux_plot)

        self.use_combustible.trace_add("write", self.combustible_callback)

        self.amb_callback()
        self.cvldlf_callback()
        self.cvldlf_consistency_callback()
        self.type_callback()
        self.ctrl_callback()
        self.material_callback()
        self.guide_callback()
        self.propellant_callback()
        self.combustible_callback()

        self.update_all()

        # self.bind("<Configure>", self.resizePlot)

        root.protocol("WM_DELETE_WINDOW", self.quit)

        self.use_theme()  # <- an update is authorized here

        self.t_lid = None

        text_handler = TextHandler(self.error_text)
        root_logger.addHandler(text_handler)
        root_logger.info("text handler attached to root logger.")

        console = logging.StreamHandler(sys.stderr)
        console.setFormatter(formatter)

        self.listener = QueueListener(self.log_queue, text_handler, console)
        self.listener.start()

        root_logger.info("text handler attached to subprocess log queue listener.")

        self.timed_loop()

    def handle_errors(self, exception: Exception, level: Literal[30, 40] = logging.ERROR):
        if self.debug.get():
            exc_type, exc_value, exc_traceback = sys.exc_info()
            root_logger.log(level, "".join(traceback.format_exception(exc_type, exc_value, exc_traceback)))
        else:
            root_logger.log(level, str(exception))

    def timed_loop(self):
        if self.process:
            self.get_value()

        if self.guide_process:
            self.get_guide()

        self.t_lid = self.root.after(100, self.timed_loop)

    def quit(self):
        if self.process:
            self.process.terminate()
            self.process.kill()

        if self.guide_process:
            self.guide_process.terminate()
            self.guide_process.kill()

        self.listener.stop()

        if self.t_lid is not None:
            self.root.after_cancel(self.t_lid)

        super().quit()

    def get_description(self):
        if self.gun is None:
            return "N/A"
        else:
            typ, cal, w = self.kwargs["typ"], self.kwargs["caliber"], self.kwargs["shot_mass"]
            return "{:} {:.0f} mm ".format(
                "{:.4g} g".format(w * 1e3) if w < 1 else "{:.4g} kg".format(w), cal * 1e3
            ) + self.get_loc_str(typ)

    def save(self):
        gun = self.gun
        if gun is None:
            messagebox.showinfo(self.get_loc_str("excTitle"), self.get_loc_str("noDataMsg"))
            return

        file_name = filedialog.asksaveasfilename(
            title=self.get_loc_str("saveLabel"),
            filetypes=(("JSON file", "*.json"),),
            defaultextension=".json",
            initialfile=filenameize(self.get_description()),
        )

        try:
            loc_val_dict = {
                loc.get_descriptive(): loc.get()
                for loc in self.locs
                if (hasattr(loc, "get_descriptive") and loc.get_descriptive())
            }

            kvs = {**loc_val_dict, "Description": self.description.get(1.0, "end").strip("\n")}
            with open(file_name, "w", encoding="utf-8") as file:
                json.dump(kvs, file, indent="\t", ensure_ascii=False, sort_keys=True)

            messagebox.showinfo(
                self.get_loc_str("sucTitle"), self.get_loc_str("savedLocMsg") + " {:}".format(file_name)
            )

        except Exception as e:
            self.handle_errors(e, logging.WARNING)

    def load_gun(self, initial_dir: str = None):
        file_name = filedialog.askopenfilename(
            title=self.get_loc_str("loadLabel"),
            filetypes=(("JSON File", "*.json"),),
            defaultextension=".json",
            initialdir=initial_dir,
        )

        if not file_name:
            return

        self.on_reset()
        self.remove_cvldlf_consistency_traces()

        try:
            loc_dict = {
                loc.get_descriptive(): loc
                for loc in self.locs
                if (hasattr(loc, "get_descriptive") and loc.get_descriptive())
            }
            with open(file_name, "r", encoding="utf-8") as file:
                file_dict = json.load(file)

            for key, value in file_dict.items():
                try:
                    if isinstance(value, bool) or isinstance(value, str):
                        loc_dict[key].set(value)
                except KeyError:
                    pass
                except ValueError:
                    loc_dict[key].reset()

            for key, value in file_dict.items():
                try:
                    if isinstance(value, float) or isinstance(value, int):
                        loc_dict[key].set(value)
                except KeyError:
                    pass
                except ValueError:
                    loc_dict[key].reset()

            if DESCRIPTION in file_dict.keys():  # update description from file.
                self.description.delete(1.0, "end")
                self.description.insert("end", file_dict[DESCRIPTION].strip("\n"))

        except Exception as e:
            self.handle_errors(e, logging.WARNING)
        finally:
            self.add_cvldlf_consistency_traces()
            self.on_calculate()

    def load_propellant(self):
        file_name = filedialog.askopenfilename(
            title=self.get_loc_str("loadLabel"), filetypes=(("Comma Separated Values File", "*.csv"),)
        )

        try:
            self.drop_prop.reset(str_obj_dict=Composition.read_file(file_name))
        except Exception as e:
            self.handle_errors(e, logging.WARNING)

    def on_reset(self):
        for loc in self.locs:
            loc.reset() if hasattr(loc, "get_descriptive") else None
        self.description.delete(1.0, "end")
        self.kwargs, self.gun, self.gun_result, self.guide = None, None, None, None
        self.update_all()

    def update_all(self):
        self.update_stats()
        self.update_table()
        self.update_geom()
        self.update_spec()
        self.update_fig_plot()
        self.update_aux_plot()
        self.update_geom_plot()
        self.update_guide_graph()

    def export_table(self):
        gun = self.gun
        if gun is None:
            messagebox.showinfo(self.get_loc_str("excTitle"), self.get_loc_str("noDataMsg"))
            return

        file_name = filedialog.asksaveasfilename(
            title=self.get_loc_str("exportLabel"),
            filetypes=(("Comma Separated File", "*.csv"),),
            defaultextension=".csv",
            initialfile=filenameize(self.get_description()),
        )
        gun_type = self.kwargs["typ"]
        column_list = self.get_loc_str("columnList")[gun_type]

        try:
            with open(file_name, "w", encoding="utf-8", newline="") as csvFile:
                csv_writer = csv.writer(csvFile, delimiter=",", quoting=csv.QUOTE_MINIMAL)
                csv_writer.writerow(column_list)
                for line in self.gun_result.get_raw_table_data():
                    csv_writer.writerow(line)

            messagebox.showinfo(
                self.get_loc_str("sucTitle"), self.get_loc_str("savedLocMsg") + " {:}".format(file_name)
            )

        except Exception as e:
            self.handle_errors(e, logging.WARNING)

    def change_lang(self):

        self.menubar.entryconfig(1, label=self.get_loc_str("dataLabel"))
        self.menubar.entryconfig(2, label=self.get_loc_str("designLabel"))
        self.menubar.entryconfig(3, label=self.get_loc_str("themeLabel"))
        self.menubar.entryconfig(4, label=self.get_loc_str("debugLabel"))

        self.design_menu.entryconfig(0, label=self.get_loc_str("saveLabel"))
        self.design_menu.entryconfig(1, label=self.get_loc_str("loadLabel"))
        self.design_menu.entryconfig(2, label=self.get_loc_str("loadPresetLabel"))
        self.design_menu.entryconfig(3, label=self.get_loc_str("resetLabel"))
        self.design_menu.entryconfig(4, label=self.get_loc_str("calcLabel"))
        self.design_menu.entryconfig(5, label=self.get_loc_str("guideLabel"))

        self.data_menu.entryconfig(0, label=self.get_loc_str("exportMain"))
        self.data_menu.entryconfig(1, label=self.get_loc_str("exportAux"))
        self.data_menu.entryconfig(2, label=self.get_loc_str("exportGeom"))
        self.data_menu.entryconfig(3, label=self.get_loc_str("exportGuide"))
        self.data_menu.entryconfig(4, label=self.get_loc_str("exportLabel"))
        self.data_menu.entryconfig(5, label=self.get_loc_str("reloadPropellant"))

        self.debug_menu.entryconfig(0, label=self.get_loc_str("enableLabel"))

        self.calc_button_tip.set(self.get_loc_str("calcButtonText"))
        # self.specs_text_tip.set(self.get_loc_str("specsText"))

        for loc_widget in self.locs:
            loc_widget.localize()

        self.calc_button.config(text=self.get_loc_str("calcLabel"))
        self.guide_button.config(text=self.get_loc_str("guideLabel"))
        self.swap_button.config(text=self.get_loc_str("swapLabel"))

        self.tab_parent.tab(self.plot_tab, text=self.get_loc_str("plotTab"))
        self.tab_parent.tab(self.table_tab, text=self.get_loc_str("tableTab"))
        self.tab_parent.tab(self.errorTab, text=self.get_loc_str("errorTab"))
        self.tab_parent.tab(self.guide_tab, text=self.get_loc_str("guideTab"))
        #
        # self.prop_tab_parent.tab(prop_frame, text=self.get_loc_str("propFrmLabel"))
        # self.prop_tab_parent.tab(self.grain_frm, text=self.get_loc_str("grainFrmLabel"))

        self.guide_plot_travel.config(text=self.get_loc_str("guidePlotTravel"))
        self.guide_plot_volume.config(text=self.get_loc_str("guidePlotVolume"))
        self.guide_plot_burnout.config(text=self.get_loc_str("guidePlotBurnout"))

        self.update_all()

        super().change_lang()

    def generate_kwargs(self):
        constrain = bool(self.solve_W_Lg.get())
        lock = bool(self.lock_Lg.get())
        optimize = bool(self.opt_lf.get())
        debug = bool(self.debug.get())
        atmosphere = bool(self.in_atmos.get())
        autofrettage = bool(self.material_is_af.get())
        # is_tc = bool(self.is_tc.get())

        gun_type = self.type_optn.get_obj()

        if self.prop is None:
            raise ValueError("Invalid propellant.")

        use_cv = self.use_cv.get_obj()

        if use_cv == USE_CV:
            chamber_volume = float(self.cv_L.get()) * 1e-3
        elif use_cv == USE_LF:
            chamber_volume = float(self.chg_kg.get()) / self.prop.rho_p / float(self.lf.get()) * 100
        elif use_cv == USE_LD:
            chamber_volume = float(self.chg_kg.get()) / float(self.ld.get())
        else:
            raise ValueError("Invalid use_cv.")

        chambrage = float(self.clr.get())
        charge_mass = float(self.chg_kg.get())
        caliber = float(self.cal_mm.get()) * 1e-3

        gun_length = float(self.tbl_mm.get()) * 1e-3
        load_fraction = float(self.lf.get()) * 1e-2

        self.kwargs = {
            "opt": optimize,
            "con": constrain,
            "deb": debug,
            "lock": lock,
            "typ": gun_type,
            "dom": self.drop_domain.get_obj(),
            "sol": self.drop_soln.get_obj(),
            "control": self.p_control.get_obj(),
            "structural_material": (
                Material(rho=float(self.material_density.get()), yield_strength=float(self.material_yield.get()) * 1e6)
                if self.use_material.get()
                else None
            ),
            "structural_safety_factor": float(self.material_ssf.get()),
            "caliber": caliber,
            "shot_mass": float(self.sht_kg.get()),
            "propellant": self.prop,
            "web": float(self.web_mm.get()) * 1e-3,
            "charge_mass": charge_mass,
            "charge_mass_ratio": float(self.chg_kg.get()) / float(self.sht_kg.get()),
            "chamber_volume": chamber_volume,
            "start_pressure": float(self.stp_MPa.get()) * 1e6,
            "length_gun": gun_length,
            "chambrage": chambrage,  # chamber expansion
            "nozzle_expansion": float(self.nozz_exp.get()),  # nozzle expansion
            "nozzle_efficiency": float(self.nozz_eff.get()) * 1e-2,  # nozzle efficiency
            "drag_coefficient": float(self.dgc.get()) * 1e-2,  # drag coefficient
            "design_pressure": float(self.p_tgt.get()) * 1e6,  # design pressure
            "design_velocity": float(self.v_tgt.get()),  # design velocity
            "tol": 10 ** -int(self.acc_exp.get()),
            "min_web": 1e-6 * float(self.min_web.get()),
            "max_length": float(self.lg_max.get()),
            "load_fraction": load_fraction,
            "step": int(self.step.get()),
            "autofrettage": autofrettage,
            "known_bore": lock,
            "min_cmr": float(self.guide_min_cmr.get()),
            "max_cmr": float(self.guide_max_cmr.get()),
            "step_cmr": float(self.guide_step_cmr.get()),
            "min_lf": float(self.guide_min_lf.get()) * 1e-2,
            "max_lf": float(self.guide_max_lf.get()) * 1e-2,
            "step_lf": float(self.guide_step_lf.get()) * 1e-2,
            "max_iter": int(self.max_iter.get()),
            "max_guess": int(self.max_guess.get()),
        }

        if atmosphere:
            self.kwargs.update(
                {
                    "ambient_p": float(self.amb_p.get()) * 1e3,
                    "ambient_rho": float(self.amb_rho.get()),
                    "ambient_gamma": float(self.amb_gamma.get()),
                }
            )
        else:
            self.kwargs.update({"ambient_p": 0, "ambient_rho": 0, "ambient_gamma": 1})

    def on_guide(self):
        self.focus()  # remove focus to force widget entry validation
        self.update()  # and wait for the event to update.
        self.generate_kwargs()

        if self.process or self.guide_process:
            return
        try:
            self.guide_process = Process(target=guide, args=(self.guide_job_queue, self.log_queue, self.kwargs))
            self.guide_process.start()
        except Exception as e:
            self.guide_process = None
            self.handle_errors(e)
            self.guide = None
            self.update_guide_graph()
        else:
            for loc in self.locs:
                try:
                    loc.inhibit()
                except AttributeError:
                    pass

            self.calc_button.config(state="disabled")
            self.guide_button.config(state="disabled")
            self.swap_button.config(state="disabled")

    def on_calculate(self):
        self.focus()  # remove focus to force widget entry validation
        self.update()  # and wait for the event to update.
        self.generate_kwargs()

        if self.process or self.guide_process:
            return
        try:
            self.process = Process(target=calculate, args=(self.job_queue, self.log_queue, self.kwargs))
            self.process.start()
        except Exception as e:
            self.process = None
            self.handle_errors(e)

            self.gun, self.gun_result = None, None
            self.update_table()
            self.update_fig_plot()
            self.update_aux_plot()
        else:
            for loc in self.locs:
                try:
                    loc.inhibit()
                except AttributeError:
                    pass

            self.calc_button.config(state="disabled")
            self.guide_button.config(state="disabled")
            self.swap_button.config(state="disabled")

    def get_value(self):
        try:
            self.kwargs, self.gun, self.gun_result = self.job_queue.get_nowait()
        except Empty:
            return

        self.process = None
        constrain = self.kwargs["con"]
        lock = self.kwargs["lock"]
        optimize = self.kwargs["opt"]

        sigfig = int(-log10(self.kwargs["tol"])) + 1
        gun = self.gun

        if gun:
            if constrain:
                webmm = round_sig(self.kwargs["web"] * 1e3, n=sigfig)
                self.web_mm.set(webmm)

                if not lock:
                    lgmm = round_sig(self.kwargs["length_gun"] * 1e3, n=sigfig)
                    self.tbl_mm.set(lgmm)

                if optimize:
                    if self.use_cv.get_obj() == USE_CV:
                        self.cv_L.set(round_sig(self.kwargs["chamber_volume"] * 1e3, n=sigfig))
                    else:
                        self.lf.set(
                            round_sig(self.kwargs["load_fraction"] * 100, n=sigfig)
                        )  # corrected "bulk" load fraction

            # p = self.tabParent.index(self.tabParent.select())  # find the currently active tab
            # if p == 2 or p == 3:  # if the guidance/logging pane is currently selected.
            #     self.tabParent.select(0)  # goto the graphing pane

        else:  # calculation results in some error
            # self.tabParent.select(3)  # goto the error pane
            pass

        self.update_stats()
        self.update_table()
        self.update_fig_plot()
        self.update_aux_plot()

        self.calc_button.config(state="normal")
        self.guide_button.config(state="normal")
        self.swap_button.config(state="normal")

        for loc in self.locs:
            try:
                loc.disinhibit()
            except AttributeError:
                pass

    def update_stats(self, *_):
        self.name.set(self.get_description())  # this would always work

        for entry in (
            self.te,
            self.be,
            self.pe,
            self.va,
            self.lx,
            self.ammo,
            self.pa,
            self.gm,
            self.pp,
            self.mv,
            self.bop,
            self.sj,
        ):
            entry.reset()
        try:
            caliber = self.kwargs["caliber"]
            eta_t, eta_b, eta_p = self.gun_result.get_eff()
            self.te.set(f"{eta_t * 100:.2f} %")
            self.be.set(f"{eta_b * 100:.2f} %")
            self.pe.set(f"{eta_p * 100:.2f} %")
            self.va.set(to_si(self.gun.v_j, unit="m/s"))
            self.lx.set(
                (
                    to_si(self.gun.l_g / caliber, unit=self.get_loc_str("calLabel")),
                    to_si((self.gun.l_g + self.gun.l_c) / caliber, unit=self.get_loc_str("calLabel")),
                )
            )
            self.ammo.set(to_si(self.gun.l_c, unit="m"))
            ps = self.gun_result.read_table_data(POINT_PEAK_SHOT).shot_pressure
            self.pa.set(to_si(ps * self.gun.s / self.gun.m, unit="m/s²"))
            tube_mass = self.gun_result.tubeMass
            self.gm.set(format_mass(tube_mass))
            peak_average_entry = self.gun_result.read_table_data(POINT_PEAK_AVG)
            peak_breech_entry = self.gun_result.read_table_data(POINT_PEAK_BREECH)
            self.pp.set(
                (
                    f"{to_si(peak_average_entry.avg_pressure, unit='Pa'):}" + self.get_loc_str("mean"),
                    f"{to_si(peak_breech_entry.breech_pressure, unit='Pa'):}" + self.get_loc_str("breech"),
                )
            )
            muzzle_entry = self.gun_result.read_table_data(POINT_EXIT)
            self.mv.set(to_si(muzzle_entry.velocity, unit="m/s"))
            try:
                burnout_entry = self.gun_result.read_table_data(POINT_BURNOUT)
                self.bop.set(f"{burnout_entry.travel / self.gun.l_g * 1e2:.2f} %")
            except ValueError:
                self.bop.set(self.get_loc_str("uncontained"))

            try:
                self.sj.set(f"{to_si(self.gun.S_j, unit='m²', unit_dim=2):}")
            except AttributeError:
                self.sj.set("N/A")

        except (AttributeError, ValueError, TypeError):
            pass

    def get_guide(self):
        try:
            self.guide = self.guide_job_queue.get_nowait()
        except Empty:
            return

        self.guide_process = None
        self.update_guide_graph()
        self.calc_button.config(state="normal")
        self.guide_button.config(state="normal")
        self.swap_button.config(state="normal")

        for loc in self.locs:
            try:
                loc.disinhibit()
            except AttributeError:
                pass

    def update_fig_plot(self, *_):
        with mpl.rc_context(CONTEXT):
            self.ax.cla()
            self.ax_p.cla()
            self.ax_v.cla()

            if self.gun:
                v_tgt = self.kwargs["design_velocity"]
                gun_type = self.kwargs["typ"]
                dom = self.kwargs["dom"]

                xs, vs, vxs, pas, pss, pbs, p0s, psis, etas = [], [], [], [], [], [], [], [], []

                for entry in self.gun_result.table_data:
                    tag, (time, l, psi, v, pb, p, ps, temp, vx, p0, eta) = "", repeat(0.0, 11)

                    tag = entry.tag
                    time = entry.time
                    l = entry.travel
                    psi = entry.burnup
                    v = entry.velocity
                    pb = entry.breech_pressure

                    if gun_type == RECOILLESS:
                        vx = entry.outflow_velocity
                        p0 = entry.stag_pressure
                        eta = entry.outflow_fraction

                    p = entry.avg_pressure
                    ps = entry.shot_pressure

                    if tag == self.p_control.get():
                        x_peak = (time * 1e3) if dom == DOMAIN_TIME else l
                        # noinspection PyTypeChecker
                        self.ax_p.spines.left.set_position(("data", x_peak))

                    if dom == DOMAIN_TIME:
                        xs.append(time * 1000)
                    elif dom == DOMAIN_LEN:
                        xs.append(l)

                    vs.append(v)
                    vxs.append(vx)
                    pas.append(p / 1e6)
                    pss.append(ps / 1e6)
                    pbs.append(pb / 1e6)
                    p0s.append(p0 / 1e6)
                    psis.append(psi)
                    etas.append(eta)

                if self.plot_breech_p.get():
                    self.ax_p.plot(
                        xs,
                        pbs,
                        c="xkcd:goldenrod",
                        label=(
                            self.get_loc_str("figBreech")
                            if gun_type == CONVENTIONAL
                            else self.get_loc_str("figNozzleP")
                        ),
                    )

                if gun_type == RECOILLESS:
                    if self.plot_stag_p.get():
                        self.ax_p.plot(xs, p0s, "seagreen", label=self.get_loc_str("figStagnation"))

                    if self.plot_nozzle_v.get():
                        self.ax_v.plot(xs, vxs, "royalblue", label=self.get_loc_str("figNozzleV"))

                    if self.plot_eta.get():
                        self.ax.plot(xs, etas, "crimson", label=self.get_loc_str("figOutflow"))

                if self.plot_avg_p.get():
                    self.ax_p.plot(xs, pas, "tab:green", label=self.get_loc_str("figAvgP"))

                if self.plot_base_p.get():
                    self.ax_p.plot(xs, pss, "yellowgreen", label=self.get_loc_str("figShotBase"))

                if gun_type == CONVENTIONAL or gun_type == RECOILLESS:
                    self.ax_p.axhline(
                        float(self.p_tgt.get()), c="tab:green", linestyle=":", label=self.get_loc_str("figTgtP")
                    )

                if self.plot_vel.get():
                    self.ax_v.plot(xs, vs, "tab:blue", label=self.get_loc_str("figShotVel"))
                self.ax_v.axhline(v_tgt, c="tab:blue", linestyle=":", label=self.get_loc_str("figTgtV"))

                if self.plot_burnup.get():
                    self.ax.plot(xs, psis, c="tab:red", label=self.get_loc_str("figPsi"))

                lines_labeled = []
                for lines, xvals in zip(
                    (self.ax_p.get_lines(), self.ax.get_lines(), self.ax_v.get_lines()),
                    (
                        (0.2 * xs[-1] + 0.8 * x_peak, xs[-1]),
                        (0, xs[-1]),
                        (x_peak, 0.2 * xs[-1] + 0.8 * x_peak),
                        (0, x_peak),
                    ),
                ):
                    labelLines(lines, align=True, xvals=xvals, outline_width=4)
                    lines_labeled.append(lines)

                self.ax.set_xlim(left=0, right=xs[-1])
                # noinspection SpellCheckingInspection
                pmax = max(pas + pbs + pss + p0s)
                self.ax_p.set(ylim=(0, pmax * 1.1))
                self.ax_v.set(ylim=(0, max(vs + vxs) * 1.15))
                self.ax.set_ylim(bottom=0, top=1.05)

                self.ax_p.yaxis.set_ticks([v for v in self.ax_p.get_yticks() if v <= pmax][1:])

                self.ax.yaxis.tick_right()
                self.ax_p.yaxis.tick_left()
                self.ax_v.yaxis.tick_left()

                self.ax.tick_params(axis="y", colors="tab:red")
                self.ax_v.tick_params(axis="y", colors="tab:blue")
                self.ax_p.tick_params(axis="y", colors="tab:green")
                self.ax.tick_params(axis="x")

                if dom == DOMAIN_TIME:
                    self.ax.set_xlabel(self.get_loc_str("figTimeDomain"))
                elif dom == DOMAIN_LEN:
                    self.ax.set_xlabel(self.get_loc_str("figLenDomain"))

                self.ax_p.set_ylabel("MPa")
                self.ax_p.yaxis.label.set_color("tab:green")

                self.ax_v.set_ylabel("m/s")
                self.ax_v.yaxis.label.set_color("tab:blue")
            else:
                pass

            self.plt_canvas.draw_idle()

    def update_aux_plot(self, *_):
        if self.gun is None:
            with mpl.rc_context(CONTEXT):
                self.aux_ax.cla()
                self.aux_ax_h.cla()
                self.aux_canvas.draw_idle()
            return

        with mpl.rc_context(CONTEXT):
            self.aux_ax.cla()
            self.aux_ax_h.cla()

            p_trace = self.gun_result.pressure_trace
            cmap = mpl.colormaps["afmhot" + ("_r" if THEMES[self.theme_name_var.get()] else "")]

            x_max, y_max, t_min, t_max = 0, 0, inf, 0
            for trace in p_trace:
                if not trace.temperature:
                    continue
                if trace.temperature > t_max:
                    t_max = trace.temperature
                elif trace.temperature < t_min:
                    t_min = trace.temperature

            for pressureTraceEntry in p_trace[::-1]:
                tag, t, trace = (
                    pressureTraceEntry.tag,
                    pressureTraceEntry.temperature,
                    pressureTraceEntry.pressure_trace,
                )

                if t:
                    v = (t - t_min) / (t_max - t_min)
                    color = cmap(v)
                else:
                    color = cmap(0.5)
                linestyle = None
                alpha = None

                x, y = zip(*[(ppp.x, ppp.p) for ppp in trace])
                y = [v * 1e-6 for v in y]
                x_max = max(x_max, max(x))
                y_max = max(y_max, max(y))

                if self.trace_press.get():
                    self.aux_ax.plot(x, y, c=color, alpha=alpha, ls=linestyle)

            self.aux_ax.set_xlim(left=0, right=x_max)
            self.aux_ax.set_ylim(bottom=0, top=y_max * 1.15)

            self.aux_ax.tick_params(axis="y", colors="tab:green")
            self.aux_ax_h.tick_params(axis="y", colors="tab:blue")
            self.aux_ax.tick_params(axis="x")

            self.aux_ax.set_xlabel(self.get_loc_str("figAuxDomain"))

            self.aux_ax.yaxis.label.set_color("tab:green")
            self.aux_ax.set_ylabel("MPa")

            self.aux_ax_h.yaxis.set_ticks_position("right")
            self.aux_ax_h.yaxis.set_label_position("right")

            self.aux_ax_h.yaxis.label.set_color("tab:blue")
            self.aux_ax_h.set_ylabel("mm")

            h_trace = self.gun_result.outline

            if h_trace is not None and self.trace_hull.get():

                x_hull, r_in, r_out = zip(*[(entry.x, entry.r_in, entry.r_ex) for entry in h_trace])
                r_in, r_out = [r * 1e3 for r in r_in], [r * 1e3 for r in r_out]

                self.aux_ax_h.plot(x_hull, r_in, c="tab:blue")
                self.aux_ax_h.plot(x_hull, r_out, c="tab:blue")

                self.aux_ax_h.fill_between(
                    x_hull, r_in, r_out, alpha=0.5 if self.trace_press.get() else 0.8, color="tab:blue"
                )

                self.aux_ax.set_xlim(left=min(x_hull))

            self.aux_ax_h.set_ylim(bottom=0)
            self.aux_canvas.draw_idle()

    def update_spec(self, *_):
        self.specs.config(state="normal")
        compo = self.drop_prop.get_obj()
        self.specs.delete("1.0", "end")

        if compo.T_v:
            self.specs.insert(
                "end",
                "{:}: {:>4.0f} K {:}\n".format(self.get_loc_str("TvDesc"), compo.T_v, self.get_loc_str("isochorDesc")),
            )
        self.specs.insert("end", "{:}: {:>4.0f} kg/m³\n".format(self.get_loc_str("densityDesc"), compo.rho_p))
        self.specs.insert("end", "{:}: {:>4.0f} kJ/kg\n".format(self.get_loc_str("force"), compo.f / 1e3))
        isp = compo.get_isp()
        self.specs.insert(
            "end", "{:}: {:>4.0f} m/s {:>3.0f} s\n".format(self.get_loc_str("vacISPDesc"), isp, isp / 9.805)
        )
        isp = compo.get_isp(50)
        self.specs.insert(
            "end",
            "{:}: {:>4.0f} m/s {:>3.0f} s\n{:}\n".format(
                self.get_loc_str("atmISPDesc"), isp, isp / 9.805, self.get_loc_str("pRatioDesc")
            ),
        )
        self.specs.insert("end", "{:}:\n".format(self.get_loc_str("brDesc")))
        for p in (100e6, 200e6, 300e6):
            self.specs.insert(
                "end",
                "{:>12}".format(to_si(compo.get_lbr(p), unit="m/s", dec=3))
                + " @ {:>12}\n".format(to_si(p, unit="Pa", dec=3)),
            )

        self.specs.insert("end", compo.desc)
        self.specs.config(state="disabled")

    def update_geom(self, *_):
        for geom, r1, r2 in zip(
            (self.main_geom.get_obj(), self.aux_geom.get_obj()),
            (self.grain_r1, self.aux_grain_r1),
            (self.grain_r2, self.aux_grain_r2),
        ):
            if geom == SimpleGeometry.SPHERE:
                r1.remove()
                r2.remove()
            elif geom == SimpleGeometry.CYLINDER or geom == SimpleGeometry.TUBE:
                r1.remove()
                r2.restore()
            else:
                r1.restore()
                r2.restore()

        for geom, web, r1, r2 in zip(
            (self.main_geom.get_obj(), self.aux_geom.get_obj()),
            (self.web_mm, None),
            (self.grain_r1, self.aux_grain_r1),
            (self.grain_r2, self.aux_grain_r2),
        ):

            if geom == SimpleGeometry.SPHERE:
                if web:
                    web.localize("diamLabel", "diaText")

            elif geom == SimpleGeometry.STRIP:
                if web:
                    web.localize("widthLabel", "widthText")
                r1.localize("htwLabel", "heightRText")
                r2.localize("ltwLabel", "stripRText")

            elif geom == SimpleGeometry.CYLINDER:
                if web:
                    web.localize("diamLabel", "diaText")
                r2.localize("ltdLabel", "cylLRText")

            elif geom == SimpleGeometry.TUBE:
                if web:
                    web.localize("arcLabel", "arcText")
                r2.localize("ltarcLabel", "ltarcText")

            else:
                if web:
                    web.localize("arcLabel", "arcText")

                r1.localize("pdtarcLabel", "pdtarcText")
                # r2.localize("ltdLabel", "perfLRText")
                r2.localize("ltarcLabel", "ltarcText")

    def update_geom_plot(self):
        with mpl.rc_context(CONTEXT):
            n = 100
            self.geom_ax.cla()
            prop = self.prop
            if prop is not None:
                zb = prop.z_b
                xs = [i / n for i in range(n + 1)]
                ys = [prop.f_sigma_z(x) for x in xs]

                xs.append(zb)
                ys.append(prop.f_sigma_z(zb))

                xs.append(xs[-1])
                ys.append(0)

                self.geom_ax.plot(xs, ys)
                self.geom_ax.grid(which="major", color="grey", linestyle="dotted")
                self.geom_ax.minorticks_on()
                self.geom_ax.set_xlim(left=0, right=min(prop.z_b, 2))
                self.geom_ax.xaxis.set_ticks([i * 0.5 for i in range(ceil(min(prop.z_b, 2) / 0.5) + 1)])
                self.geom_ax.set_ylim(bottom=0, top=max(ys))
                self.geom_ax.yaxis.set_ticks([i * 0.25 for i in range(ceil(max(ys) / 0.25) + 1)])

            self.geom_canvas.draw_idle()

    def update_table(self):
        self.tv.delete(*self.tv.get_children())
        if not self.gun:
            return

        gun_type = self.kwargs["typ"]
        table_data, error_data, use_sn, units = [], [], (), ()

        if gun_type == CONVENTIONAL:
            table_data, error_data = (
                [
                    (
                        entry.tag,
                        entry.time,
                        entry.travel,
                        entry.burnup,
                        entry.velocity,
                        entry.breech_pressure,
                        entry.avg_pressure,
                        entry.shot_pressure,
                        entry.temperature,
                    )
                    for entry in data
                ]
                for data in (self.gun_result.table_data, self.gun_result.error_data)
            )
            use_sn = False, False, False, True, False, False, False, False, True
            units = None, "s", "m", None, "m/s", "Pa", "Pa", "Pa", "K"

        elif gun_type == RECOILLESS:
            table_data, error_data = (
                [
                    (
                        entry.tag,
                        entry.time,
                        entry.travel,
                        entry.burnup,
                        entry.outflow_fraction,
                        entry.velocity,
                        entry.outflow_velocity,
                        entry.breech_pressure,
                        entry.stag_pressure,
                        entry.avg_pressure,
                        entry.shot_pressure,
                        entry.temperature,
                    )
                    for entry in data
                ]
                for data in (self.gun_result.table_data, self.gun_result.error_data)
            )

            use_sn = False, False, False, True, True, False, False, False, False, False, False, True
            units = None, "s", "m", None, None, "m/s", "m/s", "Pa", "Pa", "Pa", "Pa", "K"

        # translate the tags.
        loc_table_data = []
        for i, line in enumerate(table_data):
            loc_table_data.append((self.get_loc_str(line[0]), *line[1:]))

        loc_table_data = dot_aligned(loc_table_data, units=units, use_sn=use_sn)
        error_data = dot_aligned(error_data, units=units, use_sn=use_sn)

        column_list = self.get_loc_str("columnList")[gun_type]
        self.tv["columns"] = column_list
        self.tv["show"] = "headings"

        self.tv.tag_configure(self.get_loc_str(POINT_PEAK_STAG), foreground="#2e8b57")
        self.tv.tag_configure(self.get_loc_str(POINT_PEAK_AVG), foreground="#2ca02c")
        self.tv.tag_configure(self.get_loc_str(POINT_PEAK_BREECH), foreground="orange")
        self.tv.tag_configure(self.get_loc_str(POINT_PEAK_SHOT), foreground="yellow green")
        self.tv.tag_configure(self.get_loc_str(POINT_BURNOUT), foreground="red")
        self.tv.tag_configure(self.get_loc_str(POINT_FRACTURE), foreground="brown")
        self.tv.tag_configure(self.get_loc_str(POINT_EXIT), foreground="steel blue")
        self.tv.tag_configure(self.get_loc_str(POINT_START), foreground="steel blue")
        self.tv.tag_configure(self.get_loc_str(COMPUTE), foreground="tan")

        self.tv.tag_configure("monospace", font=self.font)
        self.tv.tag_configure("error", font=self.font, foreground="dim gray")

        # we use a fixed width font so any char will do
        font_width, _ = self.font.measure("m"), self.font.metrics("linespace")

        win_width = self.tv.winfo_width()
        width = win_width // len(self.tv["columns"])

        for i, column in enumerate(column_list):  # foreach column
            self.tv.heading(i, text=column, anchor="e")  # let the column heading = column name
            self.tv.column(column, stretch=True, width=width, minwidth=font_width * 14, anchor="e")

        for i, (row, erow) in enumerate(zip(loc_table_data, error_data)):
            self.tv.insert("", "end", str(i + 1), values=row, tags=(row[0].strip(), "monospace"))
            self.tv.insert(
                str(i + 1), "end", str(-i - 1), values=tuple("±" + e if "." in e else e for e in erow), tags="error"
            )
            self.tv.move(str(-i - 1), str(i + 1), -1)

    def update_guide_graph(self):
        style = ttk.Style(self)
        bgc = str(style.lookup("TFrame", "background"))
        fgc = str(style.lookup("TFrame", "foreground"))

        with mpl.rc_context(CONTEXT):
            self.guide_ax.cla()

            if self.guide and any(line for line in self.guide):
                load_densities = list(list(value[0] for value in line) for line in self.guide)
                charge_masses = list(list(value[1] for value in line) for line in self.guide)

                self.guide_ax.set_xlabel(self.get_loc_str("guideLDDomain"))
                self.guide_ax.set_ylabel(self.get_loc_str("guideCMDomain"))

                index = int(self.guide_index.get())
                if index in (3, 4, 5):
                    values = list(
                        list(value[index] * (1e3 if index == 4 else 1) if value[index] else inf for value in line)
                        for line in self.guide
                    )

                    if index != 5:
                        min_val = min(min(line) for line in values)
                        max_val = max(
                            max(value[index] if value[index] else 0 for value in line) for line in self.guide
                        ) * (1e3 if index == 4 else 1)
                        max_val = min(2 * min_val, max_val)
                    else:
                        min_val, max_val = 0, 1

                    self.guide_ax.pcolormesh(
                        load_densities,
                        charge_masses,
                        values,
                        shading="nearest",
                        cmap="afmhot" + ("" if THEMES[self.theme_name_var.get()] else "_r"),
                        vmin=min_val,
                        vmax=max_val,
                    )
                    threshold = 0.5 * (min_val + max_val)
                    for line in self.guide:
                        for value in line:
                            load_density, charge_mass = value[0], value[1]
                            entry = value[index]
                            entry = entry * (1e3 if index == 4 else 1) if entry else None
                            self.guide_ax.text(
                                load_density,
                                charge_mass,
                                (f"{entry:.3g}" if index != 5 else f"{entry:.3f}") if entry else "N/A",
                                color=bgc if (entry and entry < threshold) else fgc,
                                horizontalalignment="center",
                                verticalalignment="center",
                            )
                    if index == 3:
                        self.guide_ax.set_title(self.get_loc_str("guideTravelTitle"))
                    elif index == 4:
                        self.guide_ax.set_title(self.get_loc_str("guideBVTitle"))
                    else:
                        self.guide_ax.set_title(self.get_loc_str("guideBurnoutTitle"))

                else:
                    raise ValueError("unknown guidance diagram plotting control.")

                # self.guideCursor = Cursor(self.guideAx, useblit=True, color=fgc, linewidth=1)
            self.guide_canvas.draw_idle()

    def propellant_callback(self, *_):
        """
        updates the propellant object on write to the ratio entry fields
        and, on changing the propellant or geometrical specification.
        """
        geom = self.main_geom.get_obj()
        aux_geom = self.aux_geom.get_obj()
        compo = self.drop_prop.get_obj()

        try:
            self.prop = Propellant(
                composition=compo,
                main_geom=geom,
                main_r1=self.grain_r1.get(),
                main_r2=self.grain_r2.get(),
                aux_geom=aux_geom,
                web_ratio=self.aux_web_ratio.get(),
                mass_ratio=self.aux_mass_ratio.get() if self.use_aux_grain.get() else 0.0,
                aux_r1=self.aux_grain_r1.get(),
                aux_r2=self.aux_grain_r2.get(),
                combustible_force=float(self.combustible_force_kJ__kg.get() * 1e3),
                combustible_fraction=float(
                    self.combustible_mass_kg.get() / self.chg_kg.get() if self.use_combustible.get() else 0.0
                ),
                force_fudge=float(self.force_fudge.get()) * 1e-2,
            )

        except Exception as e:
            self.handle_errors(e, level=30)
            self.prop = None

        self.update_geom_plot()

    def type_callback(self, *_):

        if self.type_optn.get_obj() == CONVENTIONAL:
            self.drop_soln.enable()

            self.nozz_exp.remove()
            self.nozz_eff.remove()
            self.plot_nozzle_v.remove()

            self.plot_breech_p.localize("plotBreechP")
            self.plot_stag_p.remove()
            self.plot_eta.remove()
            self.p_control.reset({point: point for point in (POINT_PEAK_AVG, POINT_PEAK_SHOT, POINT_PEAK_BREECH)})

        else:
            self.drop_soln.set_by_obj(SOL_LAGRANGE)
            self.drop_soln.disable()

            self.nozz_exp.restore()
            self.nozz_eff.restore()
            self.plot_nozzle_v.restore()

            self.plot_breech_p.localize("plotNozzleP")
            self.plot_stag_p.restore()
            self.plot_stag_p.localize("plotStagP")
            self.plot_eta.restore()
            self.plot_eta.localize("plotEtaEsc")
            self.p_control.reset(
                {point: point for point in (POINT_PEAK_AVG, POINT_PEAK_SHOT, POINT_PEAK_STAG, POINT_PEAK_BREECH)}
            )

    def ctrl_callback(self, *_):
        if not self.solve_W_Lg.get():
            self.v_tgt.disable()
            self.p_tgt.disable()

            self.opt_lf.disable()
            self.lock_Lg.disable()
            self.min_web.disable()
            self.lg_max.disable()
            self.p_control.disable()
            self.tbl_mm.enable()

        else:
            if self.lock_Lg.get():
                self.v_tgt.disable()
                self.opt_lf.disable()
                self.tbl_mm.enable()
            else:
                self.v_tgt.enable()
                self.opt_lf.enable()
                self.tbl_mm.disable()

            if self.opt_lf.get():
                self.lock_Lg.disable()
            else:
                self.lock_Lg.enable()

            self.p_tgt.enable()

            self.min_web.enable()
            self.lg_max.enable()
            self.p_control.enable()

        for entry in (self.aux_grain_r1, self.aux_grain_r2, self.aux_web_ratio, self.aux_mass_ratio, self.aux_geom):
            entry.disable() if self.use_aux_grain.get() == 0 else entry.enable()

    def add_cvldlf_consistency_traces(self):
        for entry in (self.chg_kg, self.cv_L, self.lf, self.ld, self.chg_kg):
            entry.var.trace_add("write", self.cvldlf_consistency_callback)

    def remove_cvldlf_consistency_traces(self):
        for entry in (self.chg_kg, self.cv_L, self.lf, self.ld, self.chg_kg):
            for trace in entry.var.trace_info():
                ops, name = trace
                if "cvldlf_consistency_callback" in name:
                    entry.var.trace_remove(*trace)

    def cvldlf_consistency_callback(self, val_name="", *_):
        try:
            compo = self.drop_prop.get_obj()
            sigfig = int(self.acc_exp.get()) + 1
            w = float(self.chg_kg.get())
            cv = float(self.cv_L.get())
            lf = float(self.lf.get())
            ld = float(self.ld.get())
            rho = compo.rho_p
            self.remove_cvldlf_consistency_traces()

            if self.use_cv.get_obj() == USE_CV:
                self.lf.set(round_sig(w / cv / rho * 1e5, n=sigfig))
                self.ld.set(round_sig(w / cv * 1e3, n=sigfig))
            elif self.use_cv.get_obj() == USE_LF:
                self.cv_L.set(round_sig(w / rho / lf * 1e5, n=sigfig))
                self.ld.set(round_sig(lf * rho * 1e-2, n=sigfig))
            elif self.use_cv.get_obj() == USE_LD:
                self.lf.set(round_sig(ld / rho * 1e2, n=sigfig))
                self.cv_L.set(round_sig(w / ld * 1e3, n=sigfig))

        except (ZeroDivisionError, ValueError):
            pass
        finally:
            self.add_cvldlf_consistency_traces()

    def material_callback(self, *_):
        for entry in (self.material_yield, self.material_density, self.material_ssf, self.material_is_af):
            entry.enable() if self.use_material.get() else entry.disable()

    def amb_callback(self, *_):
        self.amb_p.enable() if self.in_atmos.get() else self.amb_p.disable()
        self.amb_rho.enable() if self.in_atmos.get() else self.amb_rho.disable()
        self.amb_gamma.enable() if self.in_atmos.get() else self.amb_gamma.disable()

    def combustible_callback(self, *_):
        self.combustible_mass_kg.enable() if self.use_combustible.get() else self.combustible_mass_kg.disable()
        (
            self.combustible_force_kJ__kg.enable()
            if self.use_combustible.get()
            else self.combustible_force_kJ__kg.disable()
        )

    def cvldlf_callback(self, *_):
        use_cv = self.use_cv.get_obj()
        self.cv_L.disable()
        self.lf.disable()
        self.ld.disable()

        if use_cv == USE_CV:
            self.cv_L.enable()
        elif use_cv == USE_LF:
            self.lf.enable()
        elif use_cv == USE_LD:
            self.ld.enable()

    def guide_callback(self, *_):
        self.update_guide_graph()

    def use_theme(self):
        style = ttk.Style(self)
        style.theme_use(self.theme_name_var.get())

        """ensure that the treeview rows are roughly the same height
        regardless of dpi. on Windows, default is Segoe UI at 9 points
        so the default row height should be around 12"""

        style.configure("Treeview", rowheight=round(12 * (FONTSIZE / 8) * self.dpi / 72.0))
        style.configure("Treeview.Heading", font=(FONTNAME, FONTSIZE))
        style.configure("TButton", font=(FONTNAME, FONTSIZE + 2, "bold"))
        style.configure("TLabelframe.Label", font=(FONTNAME, FONTSIZE + 2, "bold"))
        style.configure("SubLabelFrame.TLabelframe.Label", font=(FONTNAME, FONTSIZE + 2))
        style.configure("TCheckbutton", font=(FONTNAME, FONTSIZE))
        style.configure("TNotebook.Tab", font=(FONTNAME, FONTSIZE + 1, "bold"))

        style = ttk.Style(self)
        bgc = str(style.lookup("TFrame", "background"))
        fgc = str(style.lookup("TFrame", "foreground"))
        # noinspection SpellCheckingInspection
        fbgc = str(style.lookup("TCombobox", "fieldbackground"))

        # some widgets also needs to be manually updated
        for widget in self.force_update_on_theme_widget:
            widget.config(background=fbgc, foreground=fgc, insertbackground=fgc)

        CONTEXT.update(
            {
                "figure.facecolor": bgc,
                "figure.edgecolor": fgc,
                "axes.edgecolor": fgc,
                "axes.facecolor": fbgc,
                "axes.labelcolor": fgc,
                "text.color": fgc,
                "xtick.color": fgc,
                "ytick.color": fgc,
            }
        )

        grays = (
            [f"gray{i}" for i in [90, 80, 70]]
            if THEMES[self.theme_name_var.get()]
            else [f"gray{i}" for i in [10, 20, 30]]
        )
        self.error_text.tag_configure("gun", background=grays[0])
        self.error_text.tag_configure("recoilless", background=grays[0])
        self.error_text.tag_configure("optimize_gun", background=grays[1])
        self.error_text.tag_configure("optimize_recoilless", background=grays[1])
        self.error_text.tag_configure("guidegraph", background=grays[2])

        try:
            for fig in (self.fig, self.geom_fig, self.aux_fig, self.guide_fig):
                fig.set_facecolor(bgc)

            for ax in (self.ax, self.ax_v, self.ax_p, self.geom_ax, self.aux_ax, self.aux_ax_h, self.guide_ax):
                ax.set_facecolor(fbgc)
                for place in ("top", "bottom", "left", "right"):
                    ax.spines[place].set_color(fgc)
                ax.tick_params(axis="x", which="both", colors=fgc, labelsize=FONTSIZE, labelfontfamily=FONTNAME)
                ax.tick_params(axis="y", which="both", colors=fgc, labelsize=FONTSIZE, labelfontfamily=FONTNAME)

            self.update_geom_plot()
            self.update_fig_plot()
            self.update_aux_plot()
            self.update_guide_graph()

        except AttributeError:
            pass

        self.update_idletasks()

    def export_graph(self, save: Literal["main", "aux", "guide", "geom"]):
        file_name = filedialog.asksaveasfilename(
            title=self.get_loc_str("exportGraphLabel"),
            filetypes=(("Portable Network Graphics", "*.png"),),
            defaultextension=".png",
            initialfile=filenameize(self.get_description() + " " + save),
        )

        try:
            if save == "main":
                fig = self.fig
            elif save == "aux":
                fig = self.aux_fig
            elif save == "guide":
                fig = self.guide_fig
            elif save == "geom":
                fig = self.geom_fig
            else:
                raise ValueError("unknown save target.")

            with mpl.rc_context(CONTEXT):
                fig.savefig(file_name, transparent=True, dpi=600)

            messagebox.showinfo(self.get_loc_str("sucTitle"), self.get_loc_str("savedLocMsg") + f" {file_name:}")

        except Exception as e:
            self.handle_errors(e, logging.WARNING)

    def swap(self):
        ## this swaps the primary and auxiliary charge.
        sigfig = int(self.acc_exp.get()) + 1
        try:
            # cache the values for swapping.
            main_geom, aux_geom = self.main_geom.get(), self.aux_geom.get()
            main_r1, main_r2 = self.grain_r1.get(), self.grain_r2.get()
            aux_r1, aux_r2 = self.aux_grain_r1.get(), self.aux_grain_r2.get()
            web_ratio, mass_ratio = self.aux_web_ratio.get(), self.aux_mass_ratio.get()

            self.grain_r1.set(aux_r1)
            self.grain_r2.set(aux_r2)
            self.aux_grain_r1.set(main_r1)
            self.aux_grain_r2.set(main_r2)

            self.aux_web_ratio.set(round_sig(1.0 / web_ratio, n=sigfig))
            self.aux_mass_ratio.set(round_sig(1.0 / mass_ratio, n=sigfig))

            self.main_geom.set(aux_geom)
            self.aux_geom.set(main_geom)

        except ValueError as e:
            self.handle_errors(e)


def calculate(job_queue, log_queue, kwargs):
    root_logger.addHandler(QueueHandler(log_queue))
    root_logger.info("calculation started.")

    gun_type = kwargs["typ"]
    constrain = kwargs["con"]
    optimize = kwargs["opt"]
    lock = kwargs["lock"]

    gun, gun_result = None, None
    try:
        if constrain:
            if gun_type == CONVENTIONAL:
                constrained = Constrained(**kwargs)
            elif gun_type == RECOILLESS:
                constrained = ConstrainedRecoilless(**kwargs)
            else:
                raise ValueError("unknown gun type")

            if optimize:
                if gun_type == CONVENTIONAL or gun_type == RECOILLESS:
                    l_f, e_1, l_g = constrained.find_min_v(**kwargs)

                else:
                    raise ValueError("unknown gun type")

                kwargs.update({"load_fraction": l_f})
                chamber_volume = kwargs["charge_mass"] / kwargs["propellant"].rho_p / kwargs["load_fraction"]
                kwargs.update({"chamber_volume": chamber_volume})
            else:
                if gun_type == CONVENTIONAL or gun_type == RECOILLESS:
                    e_1, l_g = constrained.solve(**kwargs)

                else:
                    raise ValueError("unknown gun type")

            kwargs.update({"web": 2 * e_1})

            if not lock:
                kwargs.update({"length_gun": l_g})

        if gun_type == CONVENTIONAL:
            gun = Gun(**kwargs)
        elif gun_type == RECOILLESS:
            gun = Recoilless(**kwargs)
        else:
            raise ValueError("unknown gun type")

        gun_result = gun.integrate(**kwargs)
        root_logger.info("calculation concluded.")

    except Exception as e:
        gun, gun_result = None, None
        root_logger.error("exception while calculating:")
        if kwargs["deb"]:
            exc_type, exc_value, exc_traceback = sys.exc_info()
            root_logger.error("".join(traceback.format_exception(exc_type, exc_value, exc_traceback)))
        else:
            root_logger.error(str(e))
    finally:
        job_queue.put((kwargs, gun, gun_result))


def guide(guide_job_queue, log_queue, kwargs):
    root_logger.addHandler(QueueHandler(log_queue))
    root_logger.info("guidance diagram calculation started")

    guide_results = None
    try:
        guide_results = guide_graph(**kwargs)
        root_logger.info("guidance diagram calculation concluded.")

    except Exception as e:
        guide_results = None
        root_logger.error("exception while calculating guidance diagram:")
        if kwargs["deb"]:
            exc_type, exc_value, exc_traceback = sys.exc_info()
            root_logger.error("".join(traceback.format_exception(exc_type, exc_value, exc_traceback)))
        else:
            root_logger.error(str(e))

    finally:
        guide_job_queue.put(guide_results)


def main(loc: str = None):
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

    loadfont(resolvepath("ui/sarasa-fixed-sc-regular.ttf"), True, True)
    font_manager.fontManager.addfont(resolvepath("ui/sarasa-fixed-sc-regular.ttf"))

    pibs = PIBS(loc=loc)
    pibs.mainloop()


if __name__ == "__main__":
    main()
