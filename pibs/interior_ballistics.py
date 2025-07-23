import csv
import json
import locale
import logging
import multiprocessing
import platform
import sys
import tkinter.font as tk_font
import traceback
from ctypes import windll
from logging.handlers import QueueHandler, QueueListener
from math import ceil, inf, log10, pi
from multiprocessing import Process, Queue
from queue import Empty
from tkinter import (
    Frame,
    IntVar,
    Menu,
    StringVar,
    Text,
    Tk,
    filedialog,
    messagebox,
    ttk,
)
from typing import Literal

import matplotlib as mpl
from labellines import labelLines
from matplotlib import font_manager
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

from pibs.guidegraph import guideGraph

from . import __version__
from .ballistics import (
    COMPUTE,
    CONVENTIONAL,
    DOMAIN_LEN,
    DOMAIN_TIME,
    GEOMETRIES,
    MATERIALS,
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
    Constrained,
    ConstrainedRecoilless,
    GrainComp,
    Gun,
    Propellant,
    Recoilless,
    SimpleGeometry,
)
from .localized_widget import LocLabelCheck, LocalizedFrame
from .misc import (
    STRING,
    dot_aligned,
    filenameize,
    format_int_input,
    format_mass,
    loadfont,
    resolvepath,
    round_sig,
    toSI,
    validate_ce,
    validate_nn,
    validate_pi,
    unloadfont,
)
from .tip import CreateToolTip

DESCRIPTION = "Description"

root_logger = logging.getLogger()
root_logger.setLevel(logging.INFO)
handler = logging.StreamHandler(sys.stderr)
formatter = logging.Formatter("%(asctime)s - %(name)s - %(message)s")
handler.setFormatter(formatter)
root_logger.addHandler(handler)


USE_LF = "USE_LF"
USE_CV = "USE_CV"

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


# {theme_name: is_dark}
THEMES = {"awdark": 0, "awlight": 1}


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
        self.setFormatter(logging.Formatter("%(asctime)s - %(name)s - %(message)s"))
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
        self.text.after(33, append)


class PIBS(Tk):
    def __init__(self, *args, loc, **kwargs):
        super().__init__(*args, **kwargs)

        self.option_add("*tearOff", False)

        font = tk_font.Font(family=FONTNAME, size=FONTSIZE)
        self.option_add("*Font", font)
        self.iconbitmap(resolvepath("ui/logo.ico"))

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
        frame = InteriorBallisticsFrame(
            self,
            self,
            menubar,
            dpi,
            default_lang="English" if loc != "zh_CN" else "中文",
            localization_dict=STRING,
            font=font,
        )
        frame.grid(row=0, column=0, sticky="nsew")

        self.minsize(self.winfo_width(), self.winfo_height())  # set minimum size
        # self.mainloop()

    def quit(self):
        # explicitly unload the font at the end of program.
        unloadfont(resolvepath("ui/sarasa-fixed-sc-regular.ttf"), True, True)
        super().quit()


class InteriorBallisticsFrame(LocalizedFrame):

    def __init__(self, parent, root, menubar, dpi, default_lang, localization_dict, font):
        super().__init__(parent, menubar=menubar, default_lang=default_lang, localization_dict=localization_dict)

        self.font = font

        self.job_queue = Queue()
        self.guide_job_queue = Queue()
        self.log_queue = Queue()
        self.progress_queue = Queue()

        self.process = None
        self.guide_process = None

        self.dpi = dpi
        self.parent = parent
        self.root = root
        self.force_update_on_theme_widget = []

        self.menubar = menubar

        file_menu = Menu(menubar)
        menubar.add_cascade(label=self.get_loc_str("fileLabel"), menu=file_menu, underline=0)
        theme_menu = Menu(menubar)
        menubar.add_cascade(label=self.get_loc_str("themeLabel"), menu=theme_menu, underline=0)
        debug_menu = Menu(menubar)
        menubar.add_cascade(label=self.get_loc_str("debugLabel"), menu=debug_menu, underline=0)

        self.file_menu = file_menu
        self.theme_menu = theme_menu
        self.debug_menu = debug_menu

        self.theme_name_var = StringVar(value=list(THEMES.keys())[0])

        self.debug = IntVar(value=0)

        file_menu.add_command(label=self.get_loc_str("saveLabel"), command=self.save, underline=0)
        file_menu.add_command(label=self.get_loc_str("loadLabel"), command=self.load, underline=0)
        file_menu.add_separator()

        file_menu.add_command(label=self.get_loc_str("exportMain"), command=lambda: self.exportGraph(save="main"))
        file_menu.add_command(label=self.get_loc_str("exportAux"), command=lambda: self.exportGraph(save="aux"))
        file_menu.add_command(label=self.get_loc_str("exportGeom"), command=lambda: self.exportGraph(save="geom"))
        file_menu.add_command(label=self.get_loc_str("exportGuide"), command=lambda: self.exportGraph(save="guide"))

        file_menu.add_separator()
        file_menu.add_command(label=self.get_loc_str("exportLabel"), command=self.export_table)

        for themeName in THEMES.keys():
            theme_menu.add_radiobutton(
                label=themeName, variable=self.theme_name_var, value=themeName, command=self.useTheme
            )

        debug_menu.add_checkbutton(label=self.get_loc_str("enableLabel"), variable=self.debug, onvalue=1, offvalue=0)

        self.prop, self.gun, self.guide = None, None, None

        self.columnconfigure(1, weight=1)
        self.rowconfigure(1, weight=1)

        self.add_name_plate()
        self.addLeftFrm()
        self.addCenterFrm()
        self.add_right_frm()
        self.addSpecFrm()
        self.addGeomPlot()
        self.addGuidePlot()
        self.ambCallback()
        self.cvlfCallback()
        self.typeCallback()
        self.ctrlCallback()

        self.updateTable()
        self.updateSpec()
        self.updateGeom()

        # self.bind("<Configure>", self.resizePlot)

        parent.protocol("WM_DELETE_WINDOW", self.quit)

        self.useTheme()  # <- an update is authorized here

        self.tLid = None

        text_handler = TextHandler(self.error_text)
        root_logger.addHandler(text_handler)
        root_logger.info("text handler attached to root logger.")

        console = logging.StreamHandler(sys.stderr)
        console.setFormatter(logging.Formatter("%(asctime)s - %(name)s - %(message)s"))
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
        # polling function for the calculation subprocess
        try:
            p = None
            while not self.progress_queue.empty():
                pg = self.progress_queue.get_nowait()
                p = pg if p is None else max(p, pg)
            if p is not None:
                self.progress.set(p)
        except Empty:
            pass

        if self.process:
            self.getValue()

        if self.guide_process:
            self.getGuide()

        # noinspection PyTypeChecker
        self.tLid = self.after(33, self.timed_loop)

    def quit(self):
        if self.process:
            self.process.terminate()
            self.process.kill()

        if self.guide_process:
            self.guide_process.terminate()
            self.guide_process.kill()

        self.listener.stop()

        if self.tLid is not None:
            self.after_cancel(self.tLid)

        super().quit()

    def get_description(self):
        if self.gun is None:
            return "N/A"
        else:
            typ, cal, w = self.kwargs["typ"], self.kwargs["caliber"], self.kwargs["shotMass"]
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

        if file_name == "":
            messagebox.showinfo(self.get_loc_str("excTitle"), self.get_loc_str("cancelMsg"))
            return

        try:
            loc_val_dict = {
                loc.get_description(): str(loc.get()) for loc in self.locs if hasattr(loc, "getDescriptive")
            }
            kvs = {**loc_val_dict, "Description": self.description.get(1.0, "end").strip("\n")}
            with open(file_name, "w", encoding="utf-8") as file:
                json.dump(kvs, file, indent="\t", ensure_ascii=False, sort_keys=True)

            messagebox.showinfo(
                self.get_loc_str("sucTitle"), self.get_loc_str("savedLocMsg") + " {:}".format(file_name)
            )

        except Exception as e:
            messagebox.showinfo(self.get_loc_str("excTitle"), str(e))

    def load(self):
        file_name = filedialog.askopenfilename(
            title=self.get_loc_str("loadLabel"), filetypes=(("JSON File", "*.json"),), defaultextension=".json"
        )
        if file_name == "":
            messagebox.showinfo(self.get_loc_str("excTitle"), self.get_loc_str("cancelMsg"))
            return

        try:
            loc_dict = {loc.get_description(): loc for loc in self.locs if hasattr(loc, "getDescriptive")}
            with open(file_name, "r", encoding="utf-8") as file:
                file_dict = json.load(file)

            for key, value in file_dict.items():
                try:
                    loc_dict[key].set(value)
                except KeyError:
                    pass

            if DESCRIPTION in file_dict.keys():  # update description from file.
                self.description.delete(1.0, "end")
                self.description.insert("end", file_dict[DESCRIPTION].strip("\n"))

        except Exception as e:
            exc_type, exc_value, exc_traceback = sys.exc_info()

            if self.debug.get() == 1:
                err_msg = "".join(traceback.format_exception(exc_type, exc_value, exc_traceback))
            else:
                err_msg = str(e)
            messagebox.showinfo(self.get_loc_str("excTitle"), err_msg)

        self.on_calculate()

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

        if file_name == "":
            messagebox.showinfo(self.get_loc_str("excTitle"), self.get_loc_str("cancelMsg"))
            return
        try:
            with open(file_name, "w", encoding="utf-8", newline="") as csvFile:
                csv_writer = csv.writer(csvFile, delimiter=",", quoting=csv.QUOTE_MINIMAL)
                csv_writer.writerow(column_list)
                for line in self.gun_result.getRawTableData():
                    csv_writer.writerow(line)

            messagebox.showinfo(
                self.get_loc_str("sucTitle"), self.get_loc_str("savedLocMsg") + " {:}".format(file_name)
            )

        except Exception as e:
            messagebox.showinfo(self.get_loc_str("excTitle"), str(e))

    def change_lang(self):

        self.menubar.entryconfig(0, label=self.get_loc_str("fileLabel"))
        self.menubar.entryconfig(1, label=self.get_loc_str("themeLabel"))
        self.menubar.entryconfig(2, label=self.get_loc_str("debugLabel"))

        self.file_menu.entryconfig(0, label=self.get_loc_str("saveLabel"))
        self.file_menu.entryconfig(1, label=self.get_loc_str("loadLabel"))

        self.file_menu.entryconfig(3, label=self.get_loc_str("exportMain"))
        self.file_menu.entryconfig(4, label=self.get_loc_str("exportAux"))
        self.file_menu.entryconfig(5, label=self.get_loc_str("exportGeom"))
        self.file_menu.entryconfig(6, label=self.get_loc_str("exportGuide"))

        self.file_menu.entryconfig(8, label=self.get_loc_str("exportLabel"))

        self.debug_menu.entryconfig(0, label=self.get_loc_str("enableLabel"))

        self.calc_button_tip.set(self.get_loc_str("calcButtonText"))

        for locWidget in self.locs:
            locWidget.localize()

        self.calc_button.config(text=self.get_loc_str("calcLabel"))
        self.guide_button.config(text=self.get_loc_str("guideLabel"))

        self.tabParent.tab(self.plotTab, text=self.get_loc_str("plotTab"))
        self.tabParent.tab(self.tableTab, text=self.get_loc_str("tableTab"))
        self.tabParent.tab(self.errorTab, text=self.get_loc_str("errorTab"))
        self.tabParent.tab(self.guideTab, text=self.get_loc_str("guideTab"))

        self.prop_tab_parent.tab(self.prop_frm, text=self.get_loc_str("propFrmLabel"))
        self.prop_tab_parent.tab(self.grain_frm, text=self.get_loc_str("grainFrmLabel"))

        self.guidePlotTravel.config(text=self.get_loc_str("guidePlotTravel"))
        self.guidePlotVolume.config(text=self.get_loc_str("guidePlotVolume"))
        self.guidePlotBurnout.config(text=self.get_loc_str("guidePlotBurnout"))

        self.updateTable()
        self.updateGeom()
        self.updateSpec()
        self.updateFigPlot()
        self.updateAuxPlot()
        self.updateGuideGraph()

        super().change_lang()

    def add_name_plate(self):
        name_frm = self.add_localized_label_frame(self, label_loc_key="nameFrm")
        name_frm.grid(row=0, column=0, sticky="nsew", columnspan=2)

        name_frm.columnconfigure(0, weight=1)
        name_frm.rowconfigure(0, weight=1)

        self.name = StringVar(self)
        name_plate = ttk.Entry(
            name_frm, textvariable=self.name, state="disabled", justify="left", font=(FONTNAME, FONTSIZE + 2)
        )
        name_plate.grid(row=0, column=0, sticky="nsew", padx=2, pady=2, columnspan=2)

        desc_scroll = ttk.Scrollbar(name_frm, orient="vertical")
        desc_scroll.grid(row=1, column=1, sticky="nsew")
        self.description = Text(
            name_frm, wrap="word", height=3, width=80, yscrollcommand=desc_scroll.set, font=(FONTNAME, FONTSIZE)
        )
        self.description.grid(row=1, column=0, sticky="nsew")

        desc_scroll.config(command=self.description.yview)

        self.force_update_on_theme_widget.append(self.description)

    # noinspection SpellCheckingInspection
    def addLeftFrm(self):
        left_frm = ttk.Frame(self)
        left_frm.grid(row=1, column=0, sticky="nsew")
        left_frm.columnconfigure(0, weight=1)
        left_frm.rowconfigure(0, weight=1)

        par_frm = self.add_localized_label_frame(left_frm, label_loc_key="parFrmLabel")
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
        self.ld, i = self.add_localized_12_display(parent=par_frm, row=i, label_loc_key="ldLabel"), i + 2
        self.pa, i = self.add_localized_12_display(parent=par_frm, row=i, label_loc_key="paLabel"), i + 2
        self.gm, i = self.add_localized_12_display(parent=par_frm, row=i, label_loc_key="gmLabel"), i + 2

    def add_right_frm(self):
        right_frm = ttk.Frame(self)
        right_frm.grid(row=0, column=3, rowspan=2, sticky="nsew")
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
            self.add_localized_label_check(parent=env_frm, label_loc_key="atmosLabel", row=j, columnspan=3),
            j + 1,
        )
        self.in_atmos.trace_add("write", self.ambCallback)

        self.amb_p, j = (
            self.add_localized_3_input(
                parent=env_frm,
                row=j,
                label_loc_key="ambPresLabel",
                unit_text="kPa",
                default="101.325",
                validation=validation_nn,
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
            ),
            j + 1,
        )

        self.amb_gamma, j = (
            self.add_localized_3_input(
                parent=env_frm, row=j, label_loc_key="ambGamLabel", default="1.400", validation=validation_nn
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
                col=0,
                columnspan=3,
                default=0,
                label_loc_key="consButton",
                tooltip_loc_key="useConsText",
            ),
            i + 1,
        )
        self.solve_W_Lg.trace_add("write", self.ctrlCallback)

        cons_frm = self.add_localized_label_frame(
            op_frm, label_loc_key="consFrmLabel", style="SubLabelFrame.TLabelframe"
        )
        cons_frm.grid(row=i, column=0, columnspan=3, sticky="nsew", padx=2, pady=2)
        i += 1
        cons_frm.columnconfigure(0, weight=1)

        j = 0
        self.lock_Lg, j = (
            self.add_localized_label_check(
                parent=cons_frm,
                row=j,
                columnspan=3,
                default=0,
                label_loc_key="lockButton",
                tooltip_loc_key="lockText",
            ),
            j + 1,
        )
        self.lock_Lg.trace_add("write", self.ctrlCallback)

        self.opt_lf, j = (
            self.add_localized_label_check(
                parent=cons_frm,
                row=j,
                columnspan=3,
                default=0,
                label_loc_key="minTVButton",
                tooltip_loc_key="optLFText",
            ),
            j + 1,
        )
        self.opt_lf.trace_add("write", self.ctrlCallback)

        self.v_tgt, j = (
            self.add_localized_3_input(
                parent=cons_frm,
                row=j,
                label_loc_key="vTgtLabel",
                unit_text="m/s",
                default="1000.0",
                validation=validation_nn,
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
            ),
            j + 1,
        )
        # noinspection SpellCheckingInspection
        self.lg_max, j = (
            self.add_localized_3_input(
                parent=cons_frm,
                row=j,
                label_loc_key="maxLgLabel",
                unit_text="m",
                default="100.0",
                validation=validation_nn,
                color="red",
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
                col=0,
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
                col=0,
                label_loc_key="-log10(ε)",
                default="3",
                validation=validation_pi,
                formatter=format_int_input,
                color="red",
                tooltip_loc_key="tolText",
            ),
            i + 1,
        )

        self.calc_button = ttk.Button(op_frm, text=self.get_loc_str("calcLabel"), command=self.on_calculate)
        self.calc_button.grid(row=i, column=0, columnspan=3, sticky="nsew", padx=2, pady=2)
        # self.root.bind("<Return>", lambda *_: self.onCalculate())
        op_frm.rowconfigure(i, weight=1)
        i += 1

        self.progress = IntVar()
        ttk.Progressbar(op_frm, maximum=100, variable=self.progress).grid(
            row=i, column=0, columnspan=3, sticky="nsew", padx=2, pady=2
        )
        i += 1
        self.calc_button_tip = StringVar(value=self.get_loc_str("calcButtonText"))
        CreateToolTip(self.calc_button, self.calc_button_tip)

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
        self.guide_button = ttk.Button(ctrl_frm, text=self.get_loc_str("guideLabel"), command=self.onGuide)
        # self.root.bind("<space>", lambda *_: self.onGuide())
        self.guide_button.grid(row=6, column=0, columnspan=3, sticky="nsew", padx=2, pady=2)

    def generate_kwargs(self):
        constrain = self.solve_W_Lg.get() == 1
        lock = self.lock_Lg.get() == 1
        optimize = self.opt_lf.get() == 1
        debug = self.debug.get() == 1
        atmosphere = self.in_atmos.get() == 1
        autofrettage = self.is_af.get() == 1

        gun_type = self.type_optn.get_obj()

        if self.prop is None:
            raise ValueError("Invalid propellant.")

        if self.use_cv.get_obj() == USE_CV:
            chamber_volume = float(self.cvL.get()) * 1e-3
        else:
            chamber_volume = float(self.chg_kg.get()) / self.prop.rho_p / float(self.ldf.get()) * 100

        chambrage = float(self.clr.get())
        charge_mass = float(self.chg_kg.get())
        caliber = float(self.cal_mm.get()) * 1e-3

        gun_length = float(self.tbl_mm.get()) * 1e-3
        load_fraction = float(self.ldf.get()) * 1e-2

        self.kwargs = {
            "opt": optimize,
            "con": constrain,
            "deb": debug,
            "lock": lock,
            "typ": gun_type,
            "dom": self.drop_domain.get_obj(),
            "sol": self.drop_soln.get_obj(),
            "control": self.p_control.get_obj(),
            "structuralMaterial": self.drop_mat.get_obj().createMaterialAtTemp(self.drop_mat_temp.get_obj()),
            "structuralSafetyFactor": float(self.ssf.get()),
            "caliber": caliber,
            "shotMass": float(self.sht_kg.get()),
            "propellant": self.prop,
            "grainSize": float(self.web_mm.get()) * 1e-3,
            "chargeMass": charge_mass,
            "chargeMassRatio": float(self.chg_kg.get()) / float(self.sht_kg.get()),
            "chamberVolume": chamber_volume,
            "startPressure": float(self.stp_MPa.get()) * 1e6,
            "lengthGun": gun_length,
            "chambrage": chambrage,  # chamber expansion
            "nozzleExpansion": float(self.nozz_exp.get()),  # nozzle expansion
            "nozzleEfficiency": float(self.nozz_eff.get()) * 1e-2,  # nozzle efficiency
            "dragCoefficient": float(self.dgc.get()) * 1e-2,  # drag coefficient
            "designPressure": float(self.p_tgt.get()) * 1e6,  # design pressure
            "designVelocity": float(self.v_tgt.get()),  # design velocity
            "tol": 10 ** -int(self.acc_exp.get()),
            "minWeb": 1e-6 * float(self.min_web.get()),
            "maxLength": float(self.lg_max.get()),
            "loadFraction": load_fraction,
            "step": int(self.step.get()),
            "autofrettage": autofrettage,
            "knownBore": lock,
            "minCMR": float(self.guide_min_cmr.get()),
            "maxCMR": float(self.guide_max_cmr.get()),
            "stepCMR": float(self.guide_step_cmr.get()),
            "minLF": float(self.guide_min_lf.get()) * 1e-2,
            "maxLF": float(self.guide_max_lf.get()) * 1e-2,
            "stepLF": float(self.guide_step_lf.get()) * 1e-2,
        }

        if atmosphere:
            self.kwargs.update(
                {
                    "ambientP": float(self.amb_p.get()) * 1e3,
                    "ambientRho": float(self.amb_rho.get()),
                    "ambientGamma": float(self.amb_gamma.get()),
                }
            )
        else:
            self.kwargs.update({"ambientP": 0, "ambientRho": 0, "ambientGamma": 1})

    def onGuide(self):
        if self.process or self.guide_process:
            return

        self.focus()  # remove focus to force widget entry validation
        self.update()  # and wait for the event to update.

        self.guide_process = None
        try:
            self.generate_kwargs()
            self.guide_process = Process(
                target=guide, args=(self.guide_job_queue, self.progress_queue, self.log_queue, self.kwargs)
            )
            self.guide_process.start()

        except Exception as e:
            self.handle_errors(e)
            self.guide = None
            self.updateGuideGraph()
        else:
            for loc in self.locs:
                try:
                    loc.inhibit()
                except AttributeError:
                    pass

            self.calc_button.config(state="disabled")
            self.guide_button.config(state="disabled")

    def on_calculate(self):
        if self.process or self.guide_process:
            return
        self.focus()  # remove focus to force widget entry validation
        self.update()  # and wait for the event to update.

        self.process = None
        try:
            self.generate_kwargs()
            self.process = Process(
                target=calculate, args=(self.job_queue, self.progress_queue, self.log_queue, self.kwargs)
            )

            self.process.start()

        except Exception as e:
            self.handle_errors(e)

            self.gun, self.gun_result = None, None
            self.updateTable()
            self.updateFigPlot()
            self.updateAuxPlot()
        else:
            for loc in self.locs:
                try:
                    loc.inhibit()
                except AttributeError:
                    pass

            self.calc_button.config(state="disabled")
            self.guide_button.config(state="disabled")

    def getValue(self):
        try:
            self.kwargs, self.gun, self.gun_result = self.job_queue.get_nowait()
        except Empty:
            return

        self.process = None
        kwargs = self.kwargs

        constrain = self.kwargs["con"]
        lock = kwargs["lock"]
        optimize = kwargs["opt"]
        # gunType = kwargs["typ"]
        caliber = kwargs["caliber"]
        chambrage = kwargs["chambrage"]
        true_lf = kwargs["loadFraction"]

        bore_s = pi * 0.25 * caliber**2
        breech_s = bore_s * chambrage

        sigfig = int(-log10(kwargs["tol"])) + 1
        gun = self.gun

        if gun:
            if constrain:
                # noinspection SpellCheckingInspection
                webmm = round_sig(kwargs["grainSize"] * 1e3, n=sigfig)
                self.web_mm.set(webmm)

                if not lock:
                    # noinspection SpellCheckingInspection
                    lgmm = round_sig(kwargs["lengthGun"] * 1e3, n=sigfig)
                    self.tbl_mm.set(lgmm)

                if optimize:
                    if self.use_cv.get_obj() == USE_CV:
                        self.cvL.set(round_sig(kwargs["chamberVolume"] * 1e3, n=sigfig))
                    else:
                        self.ldf.set(
                            round_sig(kwargs["loadFraction"] * 100, n=sigfig)
                        )  # corrected "bulk" load fraction

            self.ld.set(toSI(true_lf * self.prop.rho_p, useSN=False, unit="kg/m³"))
            # true load density

            ps = self.gun_result.readTableData(POINT_PEAK_SHOT).shotPressure

            eta_t, eta_b, eta_p = self.gun_result.getEff()

            self.te.set(str(round(eta_t * 100, 2)) + " %")
            self.be.set(str(round(eta_b * 100, 2)) + " %")
            self.pe.set(str(round(eta_p * 100, 2)) + " %")

            self.va.set(toSI(gun.v_j, unit="m/s"))

            cartridge_len = kwargs["chamberVolume"] / breech_s  # is equivalent to chamber length

            self.lx.set(
                (
                    toSI(kwargs["lengthGun"] / caliber, unit="Cal"),
                    toSI((kwargs["lengthGun"] + cartridge_len / kwargs["chambrage"]) / kwargs["caliber"], unit="Cal"),
                )
            )

            self.ammo.set(toSI(cartridge_len, unit="m"))
            self.pa.set(toSI(ps * gun.S / gun.m, unit="m/s²"))
            self.name.set(self.get_description())

            tube_mass = self.gun_result.tubeMass
            self.gm.set(format_mass(tube_mass))

            peak_average_entry = self.gun_result.readTableData(POINT_PEAK_AVG)
            peak_breech_entry = self.gun_result.readTableData(POINT_PEAK_BREECH)

            self.pp.set(
                (
                    f"{toSI(peak_average_entry.avgPressure, unit='Pa'):}" + self.get_loc_str("mean"),
                    f"{toSI(peak_breech_entry.breechPressure, unit='Pa'):}" + self.get_loc_str("breech"),
                )
            )

            muzzle_entry = self.gun_result.readTableData(POINT_EXIT)

            self.mv.set(toSI(muzzle_entry.velocity, unit="m/s"))

            # p = self.tabParent.index(self.tabParent.select())  # find the currently active tab
            # if p == 2 or p == 3:  # if the guidance/logging pane is currently selected.
            #     self.tabParent.select(0)  # goto the graphing pane

        else:  # calculation results in some error
            # self.tabParent.select(3)  # goto the error pane
            pass

        self.updateTable()
        self.updateFigPlot()
        self.updateAuxPlot()
        self.calc_button.config(state="normal")
        self.guide_button.config(state="normal")

        for loc in self.locs:
            try:
                loc.disinhibit()
            except AttributeError:
                pass

    def getGuide(self):
        try:
            self.guide = self.guide_job_queue.get_nowait()
        except Empty:
            return

        self.guide_process = None
        self.updateGuideGraph()
        self.calc_button.config(state="normal")
        self.guide_button.config(state="normal")

        for loc in self.locs:
            try:
                loc.disinhibit()
            except AttributeError:
                pass

    # noinspection SpellCheckingInspection
    def addSpecFrm(self):
        spec_frm = self.add_localized_label_frame(self, label_loc_key="specFrmLabel")
        spec_frm.grid(row=0, column=2, rowspan=2, sticky="nsew")
        spec_frm.columnconfigure(0, weight=1)

        validation_nn, validation_ce = self.register(validate_nn), self.register(validate_ce)

        i = 0
        self.type_optn = self.add_localized_dropdown(
            parent=spec_frm,
            str_obj_dict={gun_type: gun_type for gun_type in (CONVENTIONAL, RECOILLESS)},
            desc_label_key="typeLabel",
        )
        self.type_optn.grid(row=i, column=0, sticky="nsew", padx=2, pady=2, columnspan=3)
        i += 1

        self.cal_mm, i = (
            self.add_localized_3_input(
                parent=spec_frm,
                row=i,
                label_loc_key="calLabel",
                unit_text="mm",
                default="50.0",
                validation=validation_nn,
            ),
            i + 1,
        )

        self.tbl_mm, i = (
            self.add_localized_3_input(
                parent=spec_frm,
                row=i,
                label_loc_key="tblLabel",
                unit_text="mm",
                default="3500.0",
                validation=validation_nn,
            ),
            i + 1,
        )

        self.sht_kg, i = (
            self.add_localized_3_input(
                parent=spec_frm,
                row=i,
                label_loc_key="shtLabel",
                unit_text="kg",
                default="2.0",
                validation=validation_nn,
            ),
            i + 1,
        )

        self.chg_kg, i = (
            self.add_localized_3_input(
                parent=spec_frm,
                row=i,
                label_loc_key="chgLabel",
                unit_text="kg",
                default="0.5",
                validation=validation_nn,
                tooltip_loc_key="chgText",
            ),
            i + 1,
        )

        self.prop_tab_parent = ttk.Notebook(spec_frm, padding=0)
        self.prop_tab_parent.grid(row=i, column=0, columnspan=3, sticky="nsew")
        spec_frm.rowconfigure(i, weight=1)

        self.prop_tab_parent.columnconfigure(0, weight=1)
        self.prop_tab_parent.rowconfigure(0, weight=1)

        self.prop_frm = self.add_localized_label_frame(
            self.prop_tab_parent,
            label_loc_key="propFrmLabel",
            style="SubLabelFrame.TLabelframe",
            tooltip_loc_key="specsText",
        )

        self.prop_frm.grid(row=0, column=0, sticky="nsew", padx=2, pady=2)
        self.prop_frm.rowconfigure(1, weight=1)
        self.prop_frm.columnconfigure(0, weight=1)

        self.drop_prop = self.add_localized_dropdown(
            parent=self.prop_frm,
            str_obj_dict=GrainComp.read_file(resolvepath("ballistics/resource/propellants.csv")),
            desc_label_key="propFrmLabel",
        )
        self.drop_prop.grid(row=0, column=0, columnspan=2, sticky="nsew", padx=2, pady=2)

        spec_scroll = ttk.Scrollbar(self.prop_frm, orient="vertical")
        spec_scroll.grid(row=1, rowspan=2, column=1, sticky="nsew")
        spec_h_scroll = ttk.Scrollbar(self.prop_frm, orient="horizontal")
        spec_h_scroll.grid(row=2, column=0, sticky="nsew")

        self.specs = Text(
            self.prop_frm,
            wrap="word",
            height=20,
            width=30,
            yscrollcommand=spec_scroll.set,
            xscrollcommand=spec_h_scroll.set,
            font=(FONTNAME, FONTSIZE),
        )
        self.specs.grid(row=1, column=0, sticky="nsew")
        spec_scroll.config(command=self.specs.yview)
        spec_h_scroll.config(command=self.specs.xview)

        self.force_update_on_theme_widget.append(self.specs)

        self.grain_frm = self.add_localized_label_frame(
            self.prop_tab_parent, label_loc_key="grainFrmLabel", style="SubLabelFrame.TLabelframe"
        )
        self.grain_frm.grid(row=i, column=0, columnspan=3, sticky="nsew", padx=2, pady=2)
        i += 1

        self.grain_frm.columnconfigure(0, weight=1)
        self.grain_frm.rowconfigure(0, weight=1)
        self.prop_tab_parent.add(self.prop_frm, text=self.get_loc_str("propFrmLabel"))
        self.prop_tab_parent.add(self.grain_frm, text=self.get_loc_str("grainFrmLabel"))

        self.prop_tab_parent.enable_traversal()

        j = 0
        self.geom_plot_frm = self.add_localized_label_frame(
            self.grain_frm, label_loc_key="σ(Z)", style="SubLabelFrame.TLabelframe", tooltip_loc_key="geomPlotText"
        )
        self.geom_plot_frm.grid(row=j, column=0, columnspan=3, sticky="nsew", padx=2, pady=2)
        j += 1

        self.geom = self.add_localized_dropdown(
            parent=self.grain_frm, str_obj_dict=GEOMETRIES, desc_label_key="Grain Geometry"
        )
        self.geom.grid(row=j, column=0, columnspan=3, sticky="nsew", padx=2, pady=2)
        j += 1

        self.web_mm, j = (
            self.add_localized_3_input(
                parent=self.grain_frm,
                row=j,
                desc_label_key="Web",
                unit_text="mm",
                default="1.0",
                validation=validation_nn,
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
            ),
            j + 1,
        )

        self.grain_r2, j = (
            self.add_localized_3_input(
                parent=self.grain_frm,
                row=j,
                desc_label_key="1/β",
                unit_text="x",
                default="2.5",
                validation=validation_nn,
            ),
            j + 1,
        )

        self.use_aux_grain, j = (
            self.add_localized_label_check(
                parent=self.grain_frm,
                row=j,
                label_loc_key="useAuxGrainLabel",
                columnspan=3,
                default=False,
            ),
            j + 1,
        )
        self.use_aux_grain.trace_add("write", self.ctrlCallback)
        self.use_aux_grain.trace_add("write", self.callback)

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
            ),
            k + 1,
        )

        self.aux_geom = self.add_localized_dropdown(
            parent=self.aux_grain_frm, str_obj_dict=GEOMETRIES, desc_label_key="Auxiliary Grain Geometry"
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
            ),
            k + 1,
        )

        self.aux_grain_r2, k = (
            self.add_localized_3_input(
                parent=self.aux_grain_frm,
                row=k,
                unit_text="x",
                default="2.5",
                desc_label_key="Auxiliary 1/β",
                validation=validation_nn,
            ),
            k + 1,
        )

        self.use_cv = self.add_localized_dropdown(
            parent=spec_frm, str_obj_dict={USE_CV: USE_CV, USE_LF: USE_LF}, desc_label_key="cvlfLabel"
        )
        self.use_cv.grid(row=i, column=0, columnspan=3, sticky="nsew", padx=2, pady=2)
        i += 1

        self.ldf, i = (
            self.add_localized_3_input(
                parent=spec_frm,
                row=i,
                label_loc_key="ldfLabel",
                unit_text="%",
                default="50.0",
                validation=validation_ce,
                tooltip_loc_key="ldfText",
            ),
            i + 1,
        )

        self.cvL, i = (
            self.add_localized_3_input(
                parent=spec_frm, row=i, label_loc_key="cvLabel", unit_text="L", default="1.0", validation=validation_nn
            ),
            i + 1,
        )

        self.use_cv.trace_add("write", self.cvlfCallback)
        for widget in (self.cvL, self.ldf, self.chg_kg, self.acc_exp):
            widget.trace_add("write", self.cvlfConsistencyCallback)

        self.clr, i = (
            self.add_localized_3_input(
                parent=spec_frm,
                row=i,
                label_loc_key="clrLabel",
                unit_text="x",
                default="1.5",
                validation=validation_nn,
                tooltip_loc_key="clrText",
            ),
            i + 1,
        )

        self.dgc, i = (
            self.add_localized_3_input(
                parent=spec_frm,
                row=i,
                label_loc_key="dgcLabel",
                unit_text="%",
                default="3.0",
                validation=validation_ce,
                tooltip_loc_key="dgcText",
            ),
            i + 1,
        )

        self.stp_MPa, i = (
            self.add_localized_3_input(
                parent=spec_frm,
                row=i,
                label_loc_key="stpLabel",
                unit_text="MPa",
                default="30.0",
                validation=validation_nn,
                tooltip_loc_key="stpText",
            ),
            i + 1,
        )

        self.nozz_exp, i = (
            self.add_localized_3_input(
                parent=spec_frm,
                row=i,
                label_loc_key="nozzExpLabel",
                unit_text="x",
                default="4.0",
                validation=validation_nn,
                tooltip_loc_key="nozzExpText",
            ),
            i + 1,
        )

        self.nozz_eff, i = (
            self.add_localized_3_input(
                parent=spec_frm,
                row=i,
                label_loc_key="nozzEffLabel",
                unit_text="%",
                default="92.0",
                validation=validation_ce,
                tooltip_loc_key="nozzEffText",
            ),
            i + 1,
        )

        mec_frm = self.add_localized_label_frame(spec_frm, label_loc_key="matFrmLabel")
        mec_frm.grid(row=i, column=0, sticky="nsew", columnspan=3)
        mec_frm.columnconfigure(0, weight=1)
        i += 1

        j = 0
        self.drop_mat = self.add_localized_dropdown(
            parent=mec_frm, str_obj_dict=MATERIALS, desc_label_key="matFrmLabel"
        )
        self.drop_mat.grid(row=j, column=0, columnspan=2, sticky="nsew", padx=2, pady=2)
        self.drop_mat.trace_add("write", lambda *args: self.drop_mat_temp.reset(self.drop_mat.get_obj().getTdict()))
        j += 1

        self.drop_mat_temp = self.add_localized_dropdown(
            parent=mec_frm, str_obj_dict=self.drop_mat.get_obj().getTdict(), desc_label_key="matFrmTempLabel"
        )
        self.drop_mat_temp.grid(row=j, column=0, columnspan=2, sticky="nsew", padx=2, pady=2)
        j += 1

        self.ssf, j = (
            self.add_localized_2_input(
                parent=mec_frm,
                row=j,
                col=0,
                label_loc_key="sffLabel",
                default="1.35",
                validation=validation_nn,
            ),
            j + 1,
        )

        self.is_af, i = (
            self.add_localized_label_check(parent=mec_frm, label_loc_key="afLabel", row=j, col=0, columnspan=2),
            i + 1,
        )
        spec_frm.rowconfigure(i, weight=5)

        self.drop_prop.trace_add("write", self.updateSpec)

        self.geom.trace_add("write", self.updateGeom)
        self.aux_geom.trace_add("write", self.updateGeom)

        self.type_optn.trace_add("write", self.typeCallback)

        for entry in (
            self.grain_r1,
            self.grain_r2,
            self.web_mm,
            self.aux_mass_ratio,
            self.aux_web_ratio,
            self.aux_grain_r1,
            self.aux_grain_r2,
            self.type_optn,
        ):
            entry.trace_add("write", self.callback)

    def addGeomPlot(self):
        with mpl.rc_context(CONTEXT):
            fig = Figure(dpi=96, layout="constrained")
            self.geomFig = fig
            self.geomAx = fig.add_subplot(111)

            self.geomCanvas = FigureCanvasTkAgg(fig, master=self.geom_plot_frm)
            self.geomCanvas.draw_idle()
            self.geomCanvas.get_tk_widget().place(relheight=1, relwidth=1)

    def addGuidePlot(self):
        plot_frm = self.add_localized_label_frame(self.guideTab, label_loc_key="guideFrmLabel")
        plot_frm.grid(row=0, column=0, sticky="nsew")

        with mpl.rc_context(CONTEXT):
            fig = Figure(dpi=96, layout="constrained")
            self.guide_fig = fig
            self.guide_ax = fig.add_subplot(111)
            self.guide_canvas = FigureCanvasTkAgg(fig, master=plot_frm)
            self.guide_canvas.draw_idle()
            self.guide_canvas.get_tk_widget().place(relheight=1, relwidth=1)

        control_frm = Frame(self.guideTab)
        control_frm.grid(row=1, column=0, sticky="nsew")

        self.guide_index = IntVar(value=3)
        self.guide_index.trace_add("write", callback=self.guideCallback)

        for i in range(3):
            control_frm.columnconfigure(i, weight=1)
        self.guidePlotTravel = ttk.Radiobutton(
            control_frm, text=self.get_loc_str("guidePlotTravel"), variable=self.guide_index, value=3
        )
        self.guidePlotTravel.grid(row=0, column=0, sticky="nsew")

        self.guidePlotVolume = ttk.Radiobutton(
            control_frm, text=self.get_loc_str("guidePlotVolume"), variable=self.guide_index, value=4
        )
        self.guidePlotVolume.grid(row=0, column=1, sticky="nsew")

        self.guidePlotBurnout = ttk.Radiobutton(
            control_frm, text=self.get_loc_str("guidePlotBurnout"), variable=self.guide_index, value=5
        )
        self.guidePlotBurnout.grid(row=0, column=2, sticky="nsew")

    def addCenterFrm(self):
        self.tabParent = ttk.Notebook(self, padding=0)
        self.tabParent.grid(row=1, column=1, sticky="nsew")

        self.tabParent.columnconfigure(0, weight=1)
        self.tabParent.rowconfigure(0, weight=1)

        self.plotTab = Frame(self.tabParent)
        self.plotTab.grid(row=0, column=0, stick="nsew")
        self.plotTab.rowconfigure(0, weight=2)
        self.plotTab.rowconfigure(1, weight=1)
        self.plotTab.columnconfigure(0, weight=1)

        self.tableTab = Frame(self.tabParent)
        self.tableTab.grid(row=0, column=0, sticky="nsew")
        self.tableTab.rowconfigure(0, weight=1)
        self.tableTab.columnconfigure(0, weight=1)

        self.guideTab = Frame(self.tabParent)
        self.guideTab.grid(row=0, column=0, sticky="nsew")
        self.guideTab.rowconfigure(0, weight=1)
        self.guideTab.columnconfigure(0, weight=1)

        self.errorTab = Frame(self.tabParent)
        self.errorTab.grid(row=0, column=0, sticky="nsew")
        self.errorTab.rowconfigure(0, weight=1)
        self.errorTab.columnconfigure(0, weight=1)

        self.tabParent.add(self.plotTab, text=self.get_loc_str("plotTab"))
        self.tabParent.add(self.tableTab, text=self.get_loc_str("tableTab"))
        self.tabParent.add(self.guideTab, text=self.get_loc_str("guideTab"))
        self.tabParent.add(self.errorTab, text=self.get_loc_str("errorTab"))

        self.tabParent.enable_traversal()

        self.addTblFrm()
        self.addErrFrm()

        self.addFigPlot()

    def addErrFrm(self):
        error_frm = self.add_localized_label_frame(self.errorTab, label_loc_key="errFrmLabel")
        error_frm.grid(row=0, column=0, columnspan=3, sticky="nsew")
        error_frm.columnconfigure(0, weight=1)
        error_frm.rowconfigure(0, weight=1)

        err_scroll = ttk.Scrollbar(error_frm, orient="vertical")
        err_scroll.grid(row=0, column=1, sticky="nsew")
        self.error_text = Text(
            error_frm,
            yscrollcommand=err_scroll.set,
            wrap="word",
            height=6,
            width=0,
            font=(FONTNAME, FONTSIZE),
            state="disabled",
        )
        self.error_text.grid(row=0, column=0, sticky="nsew")

        self.error_text.tag_configure(str(logging.DEBUG), foreground="tan")
        self.error_text.tag_configure(str(logging.WARNING), foreground="orange")
        self.error_text.tag_configure(str(logging.ERROR), foreground="orangered")
        self.error_text.tag_configure(str(logging.CRITICAL), foreground="red")

        err_scroll.config(command=self.error_text.yview)
        self.force_update_on_theme_widget.append(self.error_text)

    def addFigPlot(self):
        plot_frm = self.add_localized_label_frame(
            self.plotTab, label_loc_key="plotFrmLabel", tooltip_loc_key="plotText"
        )
        plot_frm.grid(row=0, column=0, sticky="nsew")

        for i in range(5):
            plot_frm.columnconfigure(i, weight=1)

        j = 1
        k = 0
        self.plot_avg_p = self.add_localized_label_check(parent=plot_frm, label_loc_key="plotAvgP", row=j, col=k)
        k += 1
        self.plot_base_p = self.add_localized_label_check(parent=plot_frm, label_loc_key="plotBaseP", row=j, col=k)
        k += 1
        self.plot_breech_p = self.add_localized_label_check(parent=plot_frm, label_loc_key="plotBreechP", row=j, col=k)
        k += 1
        self.plot_stag_p = self.add_localized_label_check(parent=plot_frm, label_loc_key="plotStagP", row=j, col=k)
        j += 1

        k = 0
        self.plot_vel = self.add_localized_label_check(parent=plot_frm, label_loc_key="plotVel", row=j, col=k)
        k += 1
        self.plot_nozzle_v = self.add_localized_label_check(parent=plot_frm, label_loc_key="plotNozzleV", row=j, col=k)
        k += 1
        self.plot_burnup = self.add_localized_label_check(parent=plot_frm, label_loc_key="plotBurnup", row=j, col=k)
        k += 1
        self.plot_eta = self.add_localized_label_check(parent=plot_frm, label_loc_key="plotEta", row=j, col=k)

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
            check.trace_add("write", self.updateFigPlot)

        plot_frm.columnconfigure(0, weight=1)
        plot_frm.rowconfigure(0, weight=1)

        plotPlaceFrm = Frame(plot_frm)
        plotPlaceFrm.grid(row=0, column=0, padx=2, pady=2, sticky="nsew", columnspan=5)

        with mpl.rc_context(CONTEXT):
            fig = Figure(dpi=96, layout="constrained")
            axes = fig.add_subplot(111)

            ax = axes
            axP = ax.twinx()
            axv = ax.twinx()

            ax.yaxis.tick_right()
            axv.yaxis.tick_left()

            ax.set_xlabel(" ")
            # noinspection PyTypeChecker
            axP.spines.right.set_position(("data", 0.5))
            axP.yaxis.set_ticks(axP.get_yticks()[1:-1:])

            self.ax, self.axP, self.axv, self.fig = ax, axP, axv, fig

            self.pltCanvas = FigureCanvasTkAgg(fig, master=plotPlaceFrm)
            self.pltCanvas.draw_idle()
            self.pltCanvas.get_tk_widget().place(relheight=1, relwidth=1)

        aux_frm = self.add_localized_label_frame(self.plotTab, label_loc_key="auxFrmLabel", tooltip_loc_key="auxText")
        aux_frm.grid(row=1, column=0, sticky="nsew")

        for i in range(2):
            aux_frm.columnconfigure(i, weight=1)

        j = 1
        k = 0
        self.traceHull = LocLabelCheck(
            parent=aux_frm,
            row=j,
            col=k,
            label_loc_key="traceHull",
            loc_func=self.get_loc_str,
            default=1,
            all_localized=self.locs,
        )

        k += 1
        self.tracePress = LocLabelCheck(
            parent=aux_frm, row=j, col=k, label_loc_key="tracePress", loc_func=self.get_loc_str, all_localized=self.locs
        )

        for check in (self.traceHull, self.tracePress):
            check.trace_add("write", self.updateAuxPlot)

        aux_frm.columnconfigure(0, weight=1)
        aux_frm.rowconfigure(0, weight=1)

        auxPlaceFrm = Frame(aux_frm)
        auxPlaceFrm.grid(row=0, column=0, padx=2, pady=2, sticky="nsew", columnspan=2)

        with mpl.rc_context(CONTEXT):
            fig = Figure(dpi=96, layout="constrained")
            axes = fig.add_subplot(111)

            ax = axes
            axH = ax.twinx()

            self.auxAx = ax  # auxiliary axes
            self.auxAxH = axH
            self.auxFig = fig
            self.auxAxH.yaxis.tick_right()
            self.auxCanvas = FigureCanvasTkAgg(fig, master=auxPlaceFrm)
            self.auxCanvas.draw_idle()
            self.auxCanvas.get_tk_widget().place(relwidth=1, relheight=1)

    # noinspection PyUnusedLocal
    def updateFigPlot(self, *args):
        with mpl.rc_context(CONTEXT):
            self.ax.cla()
            self.axP.cla()
            self.axv.cla()

            if self.gun:
                vTgt = self.kwargs["designVelocity"]
                gunType = self.kwargs["typ"]
                dom = self.kwargs["dom"]

                xs, vs = [], []

                Pas, Pss, Pbs, P0s = [], [], [], []
                psis, etas = [], []
                vxs = []

                if gunType == CONVENTIONAL:
                    for tag, t, l, psi, v, Pb, P, Ps, T in self.gun_result.getRawTableData():
                        if tag == POINT_PEAK_AVG:
                            if dom == DOMAIN_TIME:
                                xPeak = t * 1e3
                            elif dom == DOMAIN_LEN:
                                xPeak = l

                        if dom == DOMAIN_TIME:
                            xs.append(t * 1000)
                        elif dom == DOMAIN_LEN:
                            xs.append(l)

                        vs.append(v)
                        Pas.append(P / 1e6)
                        Pss.append(Ps / 1e6)
                        Pbs.append(Pb / 1e6)
                        psis.append(psi)

                elif gunType == RECOILLESS:
                    for tag, t, l, psi, v, vx, Px, P0, P, Ps, T, eta in self.gun_result.getRawTableData():
                        if tag == POINT_PEAK_AVG:
                            if dom == DOMAIN_TIME:
                                xPeak = t * 1e3
                            elif dom == DOMAIN_LEN:
                                xPeak = l

                        if dom == DOMAIN_TIME:
                            xs.append(t * 1000)
                        elif dom == DOMAIN_LEN:
                            xs.append(l)

                        # Fr = P * gun.S * (1 - gun.C_f * gun.S_j_bar)
                        vs.append(v)
                        vxs.append(vx)
                        Pas.append(P / 1e6)
                        Pss.append(Ps / 1e6)
                        Pbs.append(Px / 1e6)
                        P0s.append(P0 / 1e6)
                        # Frs.append(Fr / 1e6)
                        psis.append(psi)
                        etas.append(eta)

                # noinspection PyTypeChecker,PyUnreachableCode
                self.axP.spines.left.set_position(("data", xPeak))

                if self.plot_breech_p.get():
                    self.axP.plot(
                        xs,
                        Pbs,
                        c="xkcd:goldenrod",
                        label=(
                            self.get_loc_str("figBreech")
                            if gunType == CONVENTIONAL
                            else (
                                self.get_loc_str("figNozzleP")
                                if gunType == RECOILLESS
                                else self.get_loc_str("figBleedP")
                            )
                        ),
                    )

                if gunType == RECOILLESS:
                    if self.plot_stag_p.get():
                        self.axP.plot(xs, P0s, "seagreen", label=self.get_loc_str("figStagnation"))

                    if self.plot_nozzle_v.get():
                        self.axv.plot(xs, vxs, "royalblue", label=self.get_loc_str("figNozzleV"))

                    if self.plot_eta.get():
                        self.ax.plot(xs, etas, "crimson", label=self.get_loc_str("figOutflow"))

                if self.plot_avg_p.get():
                    self.axP.plot(xs, Pas, "tab:green", label=self.get_loc_str("figAvgP"))

                if self.plot_base_p.get():
                    self.axP.plot(xs, Pss, "yellowgreen", label=self.get_loc_str("figShotBase"))

                if gunType == CONVENTIONAL or gunType == RECOILLESS:
                    self.axP.axhline(
                        float(self.p_tgt.get()), c="tab:green", linestyle=":", label=self.get_loc_str("figTgtP")
                    )

                if self.plot_vel.get():
                    self.axv.plot(xs, vs, "tab:blue", label=self.get_loc_str("figShotVel"))
                self.axv.axhline(vTgt, c="tab:blue", linestyle=":", label=self.get_loc_str("figTgtV"))

                if self.plot_burnup.get():
                    self.ax.plot(xs, psis, c="tab:red", label=self.get_loc_str("figPsi"))

                linesLabeled = []
                for lines, xvals in zip(
                    (self.axP.get_lines(), self.ax.get_lines(), self.axv.get_lines()),
                    (
                        (0.2 * xs[-1] + 0.8 * xPeak, xs[-1]),
                        (0, xs[-1]),
                        (xPeak, 0.2 * xs[-1] + 0.8 * xPeak),
                        (0, xPeak),
                    ),
                ):
                    labelLines(lines, align=True, xvals=xvals, outline_width=4)
                    linesLabeled.append(lines)

                self.ax.set_xlim(left=0, right=xs[-1])
                # noinspection SpellCheckingInspection
                pmax = max(Pas + Pbs + Pss + P0s)
                self.axP.set(ylim=(0, pmax * 1.1))
                self.axv.set(ylim=(0, max(vs + vxs) * 1.15))
                self.ax.set_ylim(bottom=0, top=1.05)

                self.axP.yaxis.set_ticks([v for v in self.axP.get_yticks() if v <= pmax][1:])

                self.ax.yaxis.tick_right()
                self.axP.yaxis.tick_left()
                self.axv.yaxis.tick_left()

                self.ax.tick_params(axis="y", colors="tab:red")
                self.axv.tick_params(axis="y", colors="tab:blue")
                self.axP.tick_params(axis="y", colors="tab:green")
                self.ax.tick_params(axis="x")

                if dom == DOMAIN_TIME:
                    self.ax.set_xlabel(self.get_loc_str("figTimeDomain"))
                elif dom == DOMAIN_LEN:
                    self.ax.set_xlabel(self.get_loc_str("figLenDomain"))

                self.axP.set_ylabel("MPa")
                self.axP.yaxis.label.set_color("tab:green")

                self.axv.set_ylabel("m/s")
                self.axv.yaxis.label.set_color("tab:blue")
            else:
                pass

            self.pltCanvas.draw_idle()

    # noinspection PyUnusedLocal
    def updateAuxPlot(self, *args):
        if self.gun is None:
            with mpl.rc_context(CONTEXT):
                self.auxAx.cla()
                self.auxAxH.cla()
                self.auxCanvas.draw_idle()
            return

        with mpl.rc_context(CONTEXT):
            self.auxAx.cla()
            self.auxAxH.cla()

            pTrace = self.gun_result.pressureTrace
            # noinspection SpellCheckingInspection
            cmap = mpl.colormaps["afmhot" + ("_r" if THEMES[self.theme_name_var.get()] else "")]

            x_max, y_max, T_min, T_max = 0, 0, inf, 0
            for trace in pTrace:
                if not trace.T:
                    continue
                if trace.T > T_max:
                    T_max = trace.T
                elif trace.T < T_min:
                    T_min = trace.T

            for pressureTraceEntry in pTrace[::-1]:
                tag, T, trace = pressureTraceEntry.tag, pressureTraceEntry.T, pressureTraceEntry.pressureTrace

                if T:
                    v = (T - T_min) / (T_max - T_min)
                    color = cmap(v)
                else:
                    color = cmap(0.5)
                linestyle = None
                alpha = None

                x, y = zip(*[(ppp.x, ppp.p) for ppp in trace])
                y = [v * 1e-6 for v in y]
                x_max = max(x_max, max(x))
                y_max = max(y_max, max(y))

                if self.tracePress.get():
                    self.auxAx.plot(x, y, c=color, alpha=alpha, ls=linestyle)

            self.auxAx.set_xlim(left=0, right=x_max)
            self.auxAx.set_ylim(bottom=0, top=y_max * 1.15)

            self.auxAx.tick_params(axis="y", colors="tab:green")
            self.auxAxH.tick_params(axis="y", colors="tab:blue")
            self.auxAx.tick_params(axis="x")

            self.auxAx.set_xlabel(self.get_loc_str("figAuxDomain"))

            self.auxAx.yaxis.label.set_color("tab:green")
            self.auxAx.set_ylabel("MPa")

            self.auxAxH.yaxis.set_ticks_position("right")
            self.auxAxH.yaxis.set_label_position("right")

            self.auxAxH.yaxis.label.set_color("tab:blue")
            self.auxAxH.set_ylabel("mm")

            HTrace = self.gun_result.outline

            if HTrace is not None and self.traceHull.get():
                xHull, rIn, rOut = zip(*[trace.getRawLine() for trace in HTrace])
                rIn, rOut = [r * 1e3 for r in rIn], [r * 1e3 for r in rOut]

                self.auxAxH.plot(xHull, rIn, c="tab:blue")
                self.auxAxH.plot(xHull, rOut, c="tab:blue")

                self.auxAxH.fill_between(
                    xHull, rIn, rOut, alpha=0.5 if self.tracePress.get() else 0.8, color="tab:blue"
                )

                self.auxAx.set_xlim(left=min(xHull))

            self.auxAxH.set_ylim(bottom=0)
            self.auxCanvas.draw_idle()

    def addTblFrm(self):
        tbl_frm = self.add_localized_label_frame(self.tableTab, label_loc_key="tblFrmLabel")
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

    # noinspection PyUnusedLocal
    def updateSpec(self, *args):
        self.specs.config(state="normal")
        compo = self.drop_prop.get_obj()
        self.specs.delete("1.0", "end")
        width = self.specs.winfo_width() // self.font.measure("m")
        # noinspection SpellCheckingInspection
        if compo.T_v:
            self.specs.insert(
                "end",
                "{:}: {:>4.0f} K {:}\n".format(self.get_loc_str("TvDesc"), compo.T_v, self.get_loc_str("isochorDesc")),
            )
        self.specs.insert("end", "{:}: {:>4.0f} kg/m³\n".format(self.get_loc_str("densityDesc"), compo.rho_p))
        self.specs.insert("end", "{:}: {:>4.0f} kJ/kg\n".format(self.get_loc_str("force"), compo.f / 1e3))
        isp = compo.getIsp()
        self.specs.insert(
            "end", "{:}: {:>4.0f} m/s {:>3.0f} s\n".format(self.get_loc_str("vacISPDesc"), isp, isp / 9.805)
        )
        isp = compo.getIsp(50)
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
                "{:>12}".format(toSI(compo.getLBR(p), unit="m/s", dec=3))
                + " @ {:>12}\n".format(toSI(p, unit="Pa", dec=3)),
            )

        for line in compo.desc.split(","):
            self.specs.insert("end", line + "\n")
        self.specs.config(state="disabled")
        # this updates the specification description

        self.callback()
        self.cvlfConsistencyCallback()  # update the chamber volume / load fraction with current data

    # noinspection PyUnusedLocal
    def updateGeom(self, *args):

        for geom, r1, r2 in zip(
            (self.geom.get_obj(), self.aux_geom.get_obj()),
            (self.grain_r1, self.aux_grain_r1),
            (self.grain_r2, self.aux_grain_r2),
        ):
            if geom == SimpleGeometry.SPHERE:
                r1.remove()
                r2.remove()
            elif geom == SimpleGeometry.CYLINDER:
                r1.remove()
                r2.restore()
            else:
                r1.restore()
                r2.restore()

        for geom, web, r1, r2 in zip(
            (self.geom.get_obj(), self.aux_geom.get_obj()),
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

                # r1.reLocalize("", "")
                r2.localize("ltdLabel", "cylLRText")

            else:
                if web:
                    web.localize("athLabel", "arcText")

                r2.localize("ltdLabel", "perfLRText")
                # noinspection SpellCheckingInspection
                r1.localize("pdtalLabel", "pDiaRText")

        self.callback()

    def updateGeomPlot(self):
        with mpl.rc_context(CONTEXT):
            N = 10
            self.geomAx.cla()
            prop = self.prop
            if prop is not None:
                Zb = prop.Z_b
                xs = [i / N for i in range(N + 1)]
                ys = [prop.f_sigma_Z(x) for x in xs]

                xs.append(Zb)
                ys.append(prop.f_sigma_Z(Zb))

                xs.append(xs[-1])
                ys.append(0)

                self.geomAx.plot(xs, ys)
                self.geomAx.grid(which="major", color="grey", linestyle="dotted")
                self.geomAx.minorticks_on()
                self.geomAx.set_xlim(left=0, right=min(prop.Z_b, 2))
                self.geomAx.xaxis.set_ticks([i * 0.5 for i in range(ceil(min(prop.Z_b, 2) / 0.5) + 1)])
                self.geomAx.set_ylim(bottom=0, top=max(ys))
                self.geomAx.yaxis.set_ticks([i * 0.25 for i in range(ceil(max(ys) / 0.25) + 1)])

            self.geomCanvas.draw_idle()

    def updateTable(self):
        self.tv.delete(*self.tv.get_children())
        if self.gun is None:
            return

        try:
            gunType = self.kwargs["typ"]
            tableData, errorData = (self.gun_result.getRawTableData(), self.gun_result.getRawErrorData())
        except AttributeError:
            gunType = self.type_optn.get_obj()
            tableData, errorData = [], []

        locTableData = []
        for i, line in enumerate(tableData):
            locTableData.append((self.get_loc_str(line[0]), *line[1:]))

        if gunType == CONVENTIONAL:
            useSN = (False, False, False, True, False, False, False, False, True)
            units = (None, "s", "m", None, "m/s", "Pa", "Pa", "Pa", "K")
        elif gunType == RECOILLESS:
            useSN = (False, False, False, True, False, False, False, False, False, False, True, True)
            units = (None, "s", "m", None, "m/s", "m/s", "Pa", "Pa", "Pa", "Pa", "K", None)
        else:
            raise ValueError("unknown gun types")

        locTableData = dot_aligned(locTableData, units=units, use_sn=useSN)
        errorData = dot_aligned(errorData, units=units, use_sn=useSN)

        columnList = self.get_loc_str("columnList")[gunType]
        self.tv["columns"] = columnList
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
        fontWidth, _ = self.font.measure("m"), self.font.metrics("linespace")

        winWidth = self.tv.winfo_width()
        width = winWidth // len(self.tv["columns"])

        for i, column in enumerate(columnList):  # foreach column
            self.tv.heading(i, text=column, anchor="e")  # let the column heading = column name
            self.tv.column(column, stretch=True, width=width, minwidth=fontWidth * 14, anchor="e")

        # noinspection SpellCheckingInspection
        for i, (row, erow) in enumerate(zip(locTableData, errorData)):
            self.tv.insert("", "end", str(i + 1), values=row, tags=(row[0].strip(), "monospace"))
            self.tv.insert(
                str(i + 1), "end", str(-i - 1), values=tuple("±" + e if "." in e else e for e in erow), tags="error"
            )
            self.tv.move(str(-i - 1), str(i + 1), -1)

    def updateGuideGraph(self):
        style = ttk.Style(self)
        bgc = str(style.lookup("TFrame", "background"))
        fgc = str(style.lookup("TFrame", "foreground"))

        with mpl.rc_context(CONTEXT):
            self.guide_ax.cla()

            if self.guide and any(line for line in self.guide):
                loadDensities = list(list(value[0] for value in line) for line in self.guide)
                chargeMasses = list(list(value[1] for value in line) for line in self.guide)

                self.guide_ax.set_xlabel(self.get_loc_str("guideLDDomain"))
                self.guide_ax.set_ylabel(self.get_loc_str("guideCMDomain"))

                index = int(self.guide_index.get())
                if index in (3, 4, 5):
                    values = list(
                        list(value[index] * (1e3 if index == 4 else 1) if value[index] else inf for value in line)
                        for line in self.guide
                    )

                    if index != 5:
                        minVal = min(min(line) for line in values)
                        maxVal = max(
                            max(value[index] if value[index] else 0 for value in line) for line in self.guide
                        ) * (1e3 if index == 4 else 1)
                        maxVal = min(2 * minVal, maxVal)
                    else:
                        minVal, maxVal = 0, 1

                    self.guide_ax.pcolormesh(
                        loadDensities,
                        chargeMasses,
                        values,
                        shading="nearest",
                        cmap="afmhot" + ("" if THEMES[self.theme_name_var.get()] else "_r"),
                        vmin=minVal,
                        vmax=maxVal,
                    )
                    threshold = 0.5 * (minVal + maxVal)
                    for line in self.guide:
                        for value in line:
                            loadDensity, chargeMass = value[0], value[1]
                            entry = value[index]
                            entry = entry * (1e3 if index == 4 else 1) if entry else None
                            self.guide_ax.text(
                                loadDensity,
                                chargeMass,
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

    # noinspection PyUnusedLocal
    def callback(self, *args):
        """
        updates the propellant object on write to the ratio entry fields
        and, on changing the propellant or geometrical specification.
        """
        geom = self.geom.get_obj()
        auxGeom = self.aux_geom.get_obj()
        compo = self.drop_prop.get_obj()

        try:
            self.prop = Propellant(
                composition=compo,
                main_geom=geom,
                main_r1=float(self.grain_r1.get()),
                main_r2=float(self.grain_r2.get()),
                aux_geom=auxGeom,
                web_ratio=float(self.aux_web_ratio.get()),
                mass_ratio=float(self.aux_mass_ratio.get()) if self.use_aux_grain.get() else 0.0,
                aux_r1=float(self.aux_grain_r1.get()),
                aux_r2=float(self.aux_grain_r2.get()),
            )

        except Exception as e:
            self.handle_errors(e, level=30)
            self.prop = None

        self.updateGeomPlot()

    # noinspection PyUnusedLocal
    def typeCallback(self, *args):

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

    # noinspection PyUnusedLocal
    def ctrlCallback(self, *args):
        if self.solve_W_Lg.get() == 0:
            self.v_tgt.disable()
            self.p_tgt.disable()

            self.opt_lf.disable()
            self.lock_Lg.disable()
            self.min_web.disable()
            self.lg_max.disable()
            self.p_control.disable()
            self.tbl_mm.enable()

        else:
            if self.lock_Lg.get() == 1:
                self.v_tgt.disable()
                self.opt_lf.disable()
                self.tbl_mm.enable()
            else:
                self.v_tgt.enable()
                self.opt_lf.enable()
                self.tbl_mm.disable()

            if self.opt_lf.get() == 1:
                self.lock_Lg.disable()
            else:
                self.lock_Lg.enable()

            self.p_tgt.enable()

            self.min_web.enable()
            self.lg_max.enable()
            self.p_control.enable()

        for entry in [self.aux_grain_r1, self.aux_grain_r2, self.aux_web_ratio, self.aux_mass_ratio, self.aux_geom]:
            if self.use_aux_grain.get() == 0:
                entry.disable()
            else:
                entry.enable()

    # noinspection PyUnusedLocal, SpellCheckingInspection
    def cvlfConsistencyCallback(self, *args):
        prop = self.drop_prop.get()
        try:
            sigfig = int(self.acc_exp.get()) + 1
            if self.use_cv.get_obj() == USE_CV:  # use Cv
                cv = float(self.cvL.get()) / 1e3
                w = float(self.chg_kg.get())
                rho = prop.rho_p
                self.ldf.set(round_sig(w / cv / rho * 100, n=sigfig))

            else:  # using load fraction
                w = float(self.chg_kg.get())
                lf = float(self.ldf.get()) / 100
                rho = prop.rho_p
                self.cvL.set(round_sig((w / rho / lf) * 1e3, n=sigfig))

        except (ZeroDivisionError, ValueError):
            return

    # noinspection PyUnusedLocal
    def ambCallback(self, *args):
        self.amb_p.enable() if self.in_atmos.get() else self.amb_p.disable()
        self.amb_rho.enable() if self.in_atmos.get() else self.amb_rho.disable()
        self.amb_gamma.enable() if self.in_atmos.get() else self.amb_gamma.disable()

    # noinspection PyUnusedLocal, SpellCheckingInspection
    def cvlfCallback(self, *args):
        useCv = self.use_cv.get_obj() == USE_CV

        self.ldf.disable() if useCv else self.ldf.enable()
        self.cvL.enable() if useCv else self.cvL.disable()

    # noinspection PyUnusedLocal,
    def guideCallback(self, *args):
        self.updateGuideGraph()

    def useTheme(self):
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
            for fig in (self.fig, self.geomFig, self.auxFig, self.guide_fig):
                fig.set_facecolor(bgc)

            for ax in (self.ax, self.axv, self.axP, self.geomAx, self.auxAx, self.auxAxH, self.guide_ax):
                ax.set_facecolor(fbgc)
                for place in ("top", "bottom", "left", "right"):
                    ax.spines[place].set_color(fgc)
                ax.tick_params(axis="x", which="both", colors=fgc, labelsize=FONTSIZE, labelfontfamily=FONTNAME)
                ax.tick_params(axis="y", which="both", colors=fgc, labelsize=FONTSIZE, labelfontfamily=FONTNAME)

            self.updateGeomPlot()
            self.updateFigPlot()
            self.updateAuxPlot()
            self.updateGuideGraph()

        except AttributeError:
            pass

        self.update_idletasks()

    # noinspection PyPep8Naming
    def exportGraph(self, save: Literal["main", "aux", "guide", "geom"]):
        fileName = filedialog.asksaveasfilename(
            title=self.get_loc_str("exportGraphLabel"),
            filetypes=(("Portable Network Graphics", "*.png"),),
            defaultextension=".png",
            initialfile=filenameize(self.get_description() + " " + save),
        )

        if fileName == "":
            messagebox.showinfo(self.get_loc_str("excTitle"), self.get_loc_str("cancelMsg"))
            return
        try:
            if save == "main":
                fig = self.fig
            elif save == "aux":
                fig = self.auxFig
            elif save == "guide":
                fig = self.guide_fig
            elif save == "geom":
                fig = self.geomFig
            else:
                raise ValueError("unknown save destination")

            with mpl.rc_context(CONTEXT):
                fig.savefig(fileName, transparent=True, dpi=600)

            messagebox.showinfo(self.get_loc_str("sucTitle"), self.get_loc_str("savedLocMsg") + f" {fileName:}")

        except Exception as e:
            messagebox.showinfo(self.get_loc_str("excTitle"), str(e))


# noinspection PyPep8Naming
def calculate(jobQueue, progressQueue, logQueue, kwargs):
    root_logger.addHandler(QueueHandler(logQueue))
    root_logger.info("calculation started.")

    gunType = kwargs["typ"]
    constrain = kwargs["con"]
    optimize = kwargs["opt"]
    lock = kwargs["lock"]

    gun, gunResult = None, None
    try:
        if constrain:
            if gunType == CONVENTIONAL:
                constrained = Constrained(**kwargs)
            elif gunType == RECOILLESS:
                constrained = ConstrainedRecoilless(**kwargs)
            else:
                raise ValueError("unknown gun type")

            if optimize:
                if gunType == CONVENTIONAL or gunType == RECOILLESS:
                    l_f, e_1, l_g = constrained.findMinV(**kwargs, progressQueue=progressQueue)

                else:
                    raise ValueError("unknown gun type")

                kwargs.update({"loadFraction": l_f})
                chamberVolume = kwargs["chargeMass"] / kwargs["propellant"].rho_p / kwargs["loadFraction"]
                kwargs.update({"chamberVolume": chamberVolume})
            else:
                if gunType == CONVENTIONAL or gunType == RECOILLESS:
                    e_1, l_g = constrained.solve(**kwargs, progressQueue=progressQueue)

                else:
                    raise ValueError("unknown gun type")

            kwargs.update({"grainSize": 2 * e_1})

            if not lock:
                kwargs.update({"lengthGun": l_g})

        if gunType == CONVENTIONAL:
            gun = Gun(**kwargs)
        elif gunType == RECOILLESS:
            gun = Recoilless(**kwargs)
        else:
            raise ValueError("unknown gun type")

        gunResult = gun.integrate(**kwargs, progressQueue=progressQueue)
        root_logger.info("calculation concluded successfully.")

    except Exception as e:
        gun, gunResult = None, None
        root_logger.error("exception while calculating:")
        if kwargs["deb"]:
            exc_type, exc_value, exc_traceback = sys.exc_info()
            root_logger.error("".join(traceback.format_exception(exc_type, exc_value, exc_traceback)))
        else:
            root_logger.error(str(e))
    finally:
        jobQueue.put((kwargs, gun, gunResult))


# noinspection PyPep8Naming
def guide(guideJobQueue, progressQueue, logQueue, kwargs):
    root_logger.addHandler(QueueHandler(logQueue))
    root_logger.info("guidance diagram calculation started")

    guideResults = None
    try:
        guideResults = guideGraph(**kwargs)
        root_logger.info("guidance diagram calculation concluded successfully.")

    except Exception as e:
        guideResults = None
        root_logger.error("exception while calculating guidance diagram:")
        if kwargs["deb"]:
            exc_type, exc_value, exc_traceback = sys.exc_info()
            root_logger.error("".join(traceback.format_exception(exc_type, exc_value, exc_traceback)))
        else:
            root_logger.error(str(e))

    finally:
        guideJobQueue.put(guideResults)


# noinspection SpellCheckingInspection
def main(loc: str = None):
    multiprocessing.freeze_support()
    root_logger.info("Initializing")

    # this allows us to set our own taskbar icon
    windll.shell32.SetCurrentProcessExplicitAppUserModelID("Phoenix.Interior.Ballistics.Solver")

    # this tells windows that our program will handle scaling ourselves
    winRelease = platform.release()
    if winRelease in ("8", "10", "11"):
        windll.shcore.SetProcessDpiAwareness(1)
    elif winRelease in ("7", "Vista"):
        windll.user32.SetProcessDPIAware()
    else:
        print("Unknown release: ", winRelease, ", skipping DPI handling")

    if not loc:
        loc = locale.windows_locale[windll.kernel32.GetUserDefaultUILanguage()]

    loadfont(resolvepath("ui/sarasa-fixed-sc-regular.ttf"), True, True)
    font_manager.fontManager.addfont(resolvepath("ui/sarasa-fixed-sc-regular.ttf"))

    pibs = PIBS(loc=loc)
    pibs.mainloop()


if __name__ == "__main__":
    main()
