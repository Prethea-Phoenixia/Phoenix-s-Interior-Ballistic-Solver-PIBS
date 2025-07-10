import csv
import json
import locale
import logging
import multiprocessing
import platform
import sys
import tkinter.font as tkFont
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
from .localized_widget import (
    Loc2Input,
    Loc3Input,
    Loc12Disp,
    Loc122Disp,
    LocDropdown,
    LocLabelCheck,
    LocLabelFrame,
)
from .misc import (
    STRING,
    dot_aligned,
    filenameize,
    formatIntInput,
    formatMass,
    loadfont,
    resolvepath,
    roundSig,
    toSI,
    validateCE,
    validateNN,
    validatePI,
)
from .tip import CreateToolTip

DESCRIPTION = "Description"

rootLogger = logging.getLogger()
rootLogger.setLevel(logging.INFO)
handler = logging.StreamHandler(sys.stderr)
formatter = logging.Formatter("%(asctime)s - %(name)s - %(message)s")
handler.setFormatter(formatter)
rootLogger.addHandler(handler)


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


# noinspection PyAttributeOutsideInit,PyPep8Naming
class InteriorBallisticsFrame(Frame):

    def __init__(self, parent, root, menubar, dpi, defaultLang=None):
        super().__init__(parent)

        self.grid(row=0, column=0, sticky="nsew")
        parent.rowconfigure(0, weight=1)
        parent.columnconfigure(0, weight=1)
        self.langVar = StringVar(value=list(STRING.keys())[0] if defaultLang is None else defaultLang)

        self.locs = []

        default_font = tkFont.Font(family=FONTNAME, size=FONTSIZE)
        self.option_add("*Font", default_font)

        self.jobQueue = Queue()
        self.guideJobQueue = Queue()
        self.logQueue = Queue()
        self.progressQueue = Queue()

        self.process = None
        self.guideProcess = None

        self.dpi = dpi
        self.parent = parent
        self.root = root
        self.forceUpdOnThemeWidget = []

        self.menubar = menubar

        fileMenu = Menu(menubar)
        menubar.add_cascade(label=self.getLocStr("fileLabel"), menu=fileMenu, underline=0)
        themeMenu = Menu(menubar)
        menubar.add_cascade(label=self.getLocStr("themeLabel"), menu=themeMenu, underline=0)
        debugMenu = Menu(menubar)
        menubar.add_cascade(label=self.getLocStr("debugLabel"), menu=debugMenu, underline=0)
        langMenu = Menu(menubar)
        menubar.add_cascade(label="Lang 语言", menu=langMenu, underline=0)

        self.fileMenu = fileMenu
        self.themeMenu = themeMenu
        self.debugMenu = debugMenu

        self.themeNameVar = StringVar(value=list(THEMES.keys())[0])

        self.debug = IntVar(value=0)

        fileMenu.add_command(label=self.getLocStr("saveLabel"), command=self.save, underline=0)
        fileMenu.add_command(label=self.getLocStr("loadLabel"), command=self.load, underline=0)
        fileMenu.add_separator()

        fileMenu.add_command(label=self.getLocStr("exportMain"), command=lambda: self.exportGraph(save="main"))
        fileMenu.add_command(label=self.getLocStr("exportAux"), command=lambda: self.exportGraph(save="aux"))
        fileMenu.add_command(label=self.getLocStr("exportGeom"), command=lambda: self.exportGraph(save="geom"))
        fileMenu.add_command(label=self.getLocStr("exportGuide"), command=lambda: self.exportGraph(save="guide"))

        fileMenu.add_separator()
        fileMenu.add_command(label=self.getLocStr("exportLabel"), command=self.exportTable)

        for themeName in THEMES.keys():
            themeMenu.add_radiobutton(label=themeName, variable=self.themeNameVar, value=themeName, command=self.useTheme)

        debugMenu.add_checkbutton(label=self.getLocStr("enableLabel"), variable=self.debug, onvalue=1, offvalue=0)

        for lang in STRING.keys():
            langMenu.add_radiobutton(label=lang, variable=self.langVar, value=lang, command=self.changeLang)

        self.prop, self.gun, self.guide = None, None, None

        self.columnconfigure(1, weight=1)
        self.rowconfigure(1, weight=1)

        self.addNamePlate()
        self.addLeftFrm()
        self.addCenterFrm()
        self.addRightFrm()
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

        textHandler = TextHandler(self.errorText)
        rootLogger.addHandler(textHandler)
        rootLogger.info("text handler attached to root logger.")

        console = logging.StreamHandler(sys.stderr)
        console.setFormatter(logging.Formatter("%(asctime)s - %(name)s - %(message)s"))
        self.listener = QueueListener(self.logQueue, textHandler, console)
        self.listener.start()

        rootLogger.info("text handler attached to subprocess log queue listener.")

        self.timedLoop()

    def timedLoop(self):
        # polling function for the calculation subprocess
        try:
            p = None
            while not self.progressQueue.empty():
                pg = self.progressQueue.get_nowait()
                p = pg if p is None else max(p, pg)
            if p is not None:
                self.progress.set(p)
        except Empty:
            pass

        if self.process:
            self.getValue()

        if self.guideProcess:
            self.getGuide()

        # noinspection PyTypeChecker
        self.tLid = self.after(33, self.timedLoop)

    def quit(self):
        if self.process:
            self.process.terminate()
            self.process.kill()

        if self.guideProcess:
            self.guideProcess.terminate()
            self.guideProcess.kill()

        self.listener.stop()

        if self.tLid is not None:
            self.after_cancel(self.tLid)

        super().quit()

    def getDescriptive(self):
        if self.gun is None:
            return "N/A"
        else:
            kwargs = self.kwargs
            typ = kwargs["typ"]
            cal = kwargs["caliber"]
            w = kwargs["shotMass"]
            return "{:} {:.0f} mm ".format(
                "{:.4g} g".format(w * 1e3) if w < 1 else "{:.4g} kg".format(w), cal * 1e3
            ) + self.getLocStr(typ)

    def save(self):
        gun = self.gun
        if gun is None:
            messagebox.showinfo(self.getLocStr("excTitle"), self.getLocStr("noDataMsg"))
            return

        fileName = filedialog.asksaveasfilename(
            title=self.getLocStr("saveLabel"),
            filetypes=(("JSON file", "*.json"),),
            defaultextension=".json",
            initialfile=filenameize(self.getDescriptive()),
        )

        if fileName == "":
            messagebox.showinfo(self.getLocStr("excTitle"), self.getLocStr("cancelMsg"))
            return

        try:
            locValDict = {loc.getDescriptive(): str(loc.get()) for loc in self.locs if hasattr(loc, "getDescriptive")}
            kvs = {**locValDict, "Description": self.description.get(1.0, "end").strip("\n")}
            with open(fileName, "w", encoding="utf-8") as file:
                json.dump(kvs, file, indent="\t", ensure_ascii=False, sort_keys=True)

            messagebox.showinfo(self.getLocStr("sucTitle"), self.getLocStr("savedLocMsg") + " {:}".format(fileName))

        except Exception as e:
            messagebox.showinfo(self.getLocStr("excTitle"), str(e))

    def load(self):
        fileName = filedialog.askopenfilename(
            title=self.getLocStr("loadLabel"), filetypes=(("JSON File", "*.json"),), defaultextension=".json"
        )
        if fileName == "":
            messagebox.showinfo(self.getLocStr("excTitle"), self.getLocStr("cancelMsg"))
            return

        try:
            locDict = {loc.getDescriptive(): loc for loc in self.locs if hasattr(loc, "getDescriptive")}
            with open(fileName, "r", encoding="utf-8") as file:
                fileDict = json.load(file)

            for key, value in fileDict.items():
                try:
                    locDict[key].set(value)
                except KeyError:
                    pass

            if DESCRIPTION in fileDict.keys():  # update description from file.
                self.description.delete(1.0, "end")
                self.description.insert("end", fileDict[DESCRIPTION].strip("\n"))

        except Exception as e:
            exc_type, exc_value, exc_traceback = sys.exc_info()

            if self.debug.get() == 1:
                errMsg = "".join(traceback.format_exception(exc_type, exc_value, exc_traceback))
            else:
                errMsg = str(e)
            messagebox.showinfo(self.getLocStr("excTitle"), errMsg)

        self.onCalculate()

    def exportTable(self):
        gun = self.gun
        if gun is None:
            messagebox.showinfo(self.getLocStr("excTitle"), self.getLocStr("noDataMsg"))
            return

        fileName = filedialog.asksaveasfilename(
            title=self.getLocStr("exportLabel"),
            filetypes=(("Comma Separated File", "*.csv"),),
            defaultextension=".csv",
            initialfile=filenameize(self.getDescriptive()),
        )
        gunType = self.kwargs["typ"]
        columnList = self.getLocStr("columnList")[gunType]

        if fileName == "":
            messagebox.showinfo(self.getLocStr("excTitle"), self.getLocStr("cancelMsg"))
            return
        try:
            with open(fileName, "w", encoding="utf-8", newline="") as csvFile:
                csvWriter = csv.writer(csvFile, delimiter=",", quoting=csv.QUOTE_MINIMAL)
                csvWriter.writerow(columnList)
                for line in self.gunResult.getRawTableData():
                    csvWriter.writerow(line)

            messagebox.showinfo(self.getLocStr("sucTitle"), self.getLocStr("savedLocMsg") + " {:}".format(fileName))

        except Exception as e:
            messagebox.showinfo(self.getLocStr("excTitle"), str(e))

    def changeLang(self):
        self.menubar.entryconfig(0, label=self.getLocStr("fileLabel"))
        self.menubar.entryconfig(1, label=self.getLocStr("themeLabel"))
        self.menubar.entryconfig(2, label=self.getLocStr("debugLabel"))

        self.fileMenu.entryconfig(0, label=self.getLocStr("saveLabel"))
        self.fileMenu.entryconfig(1, label=self.getLocStr("loadLabel"))

        self.fileMenu.entryconfig(3, label=self.getLocStr("exportMain"))
        self.fileMenu.entryconfig(4, label=self.getLocStr("exportAux"))
        self.fileMenu.entryconfig(5, label=self.getLocStr("exportGeom"))
        self.fileMenu.entryconfig(6, label=self.getLocStr("exportGuide"))

        self.fileMenu.entryconfig(8, label=self.getLocStr("exportLabel"))

        self.debugMenu.entryconfig(0, label=self.getLocStr("enableLabel"))

        self.calcButtonTip.set(self.getLocStr("calcButtonText"))

        for locWidget in self.locs:
            locWidget.reLocalize()

        self.calcButton.config(text=self.getLocStr("calcLabel"))
        self.guideButton.config(text=self.getLocStr("guideLabel"))

        self.tabParent.tab(self.plotTab, text=self.getLocStr("plotTab"))
        self.tabParent.tab(self.tableTab, text=self.getLocStr("tableTab"))
        self.tabParent.tab(self.errorTab, text=self.getLocStr("errorTab"))
        self.tabParent.tab(self.guideTab, text=self.getLocStr("guideTab"))

        self.propTabParent.tab(self.propFrm, text=self.getLocStr("propFrmLabel"))
        self.propTabParent.tab(self.grainFrm, text=self.getLocStr("grainFrmLabel"))

        self.guidePlotTravel.config(text=self.getLocStr("guidePlotTravel"))
        self.guidePlotVolume.config(text=self.getLocStr("guidePlotVolume"))
        self.guidePlotBurnout.config(text=self.getLocStr("guidePlotBurnout"))

        self.updateTable()
        self.updateGeom()
        self.updateSpec()
        self.updateFigPlot()
        self.updateAuxPlot()
        self.updateGuideGraph()

    def getLocStr(self, name, forceDefault: bool = False):
        try:

            return STRING["English" if forceDefault else self.langVar.get()][name]
        except KeyError:
            try:
                return STRING["English"][name]
            except KeyError:
                return name

    def addNamePlate(self):
        nameFrm = LocLabelFrame(self, locKey="nameFrm", locFunc=self.getLocStr, allLLF=self.locs)
        nameFrm.grid(row=0, column=0, sticky="nsew", columnspan=2)

        nameFrm.columnconfigure(0, weight=1)
        nameFrm.rowconfigure(0, weight=1)
        self.name = StringVar(self)
        namePlate = ttk.Entry(
            nameFrm, textvariable=self.name, state="disabled", justify="left", font=(FONTNAME, FONTSIZE + 2)
        )
        namePlate.grid(row=0, column=0, sticky="nsew", padx=2, pady=2, columnspan=2)

        descScroll = ttk.Scrollbar(nameFrm, orient="vertical")
        descScroll.grid(row=1, column=1, sticky="nsew")
        self.description = Text(
            nameFrm, wrap="word", height=3, width=80, yscrollcommand=descScroll.set, font=(FONTNAME, FONTSIZE)
        )
        self.description.grid(row=1, column=0, sticky="nsew")

        descScroll.config(command=self.description.yview)

        self.forceUpdOnThemeWidget.append(self.description)

    def addLeftFrm(self):
        leftFrm = ttk.Frame(self)
        leftFrm.grid(row=1, column=0, sticky="nsew")
        leftFrm.columnconfigure(0, weight=1)
        leftFrm.rowconfigure(0, weight=1)

        parFrm = LocLabelFrame(leftFrm, locKey="parFrmLabel", locFunc=self.getLocStr, allLLF=self.locs)
        parFrm.grid(row=0, column=0, sticky="nsew")
        parFrm.columnconfigure(0, weight=1)

        i = 0
        self.ammo = Loc12Disp(
            parent=parFrm,
            row=i,
            labelLocKey="ammoLabel",
            unitText="",
            locFunc=self.getLocStr,
            allDisps=self.locs,
        )
        i += 2

        self.pp = Loc122Disp(
            parent=parFrm,
            row=i,
            labelLocKey="ppLabel",
            tooltipLocKey="ppText",
            locFunc=self.getLocStr,
            allDisps=self.locs,
        )
        i += 3

        self.lx = Loc122Disp(
            parent=parFrm,
            row=i,
            labelLocKey="lxLabel",
            tooltipLocKey="calLxText",
            locFunc=self.getLocStr,
            allDisps=self.locs,
        )
        i += 3

        self.mv = Loc12Disp(parent=parFrm, row=i, labelLocKey="mvLabel", locFunc=self.getLocStr, allDisps=self.locs)
        i += 2

        self.va = Loc12Disp(
            parent=parFrm,
            row=i,
            labelLocKey="vaLabel",
            tooltipLocKey="vinfText",
            locFunc=self.getLocStr,
            allDisps=self.locs,
        )
        i += 2
        self.te = Loc12Disp(
            parent=parFrm,
            row=i,
            labelLocKey="teffLabel",
            tooltipLocKey="teffText",
            locFunc=self.getLocStr,
            allDisps=self.locs,
        )
        i += 2
        # noinspection SpellCheckingInspection
        self.be = Loc12Disp(
            parent=parFrm,
            row=i,
            labelLocKey="beffLabel",
            tooltipLocKey="beffText",
            locFunc=self.getLocStr,
            allDisps=self.locs,
        )
        i += 2
        # noinspection SpellCheckingInspection
        self.pe = Loc12Disp(
            parent=parFrm,
            row=i,
            labelLocKey="peffLabel",
            tooltipLocKey="peffText",
            locFunc=self.getLocStr,
            allDisps=self.locs,
        )
        i += 2

        self.ld = Loc12Disp(parent=parFrm, row=i, labelLocKey="ldLabel", locFunc=self.getLocStr, allDisps=self.locs)
        i += 2

        self.pa = Loc12Disp(parent=parFrm, row=i, labelLocKey="paLabel", locFunc=self.getLocStr, allDisps=self.locs)
        i += 2

        self.gm = Loc12Disp(parent=parFrm, row=i, labelLocKey="gmLabel", locFunc=self.getLocStr, allDisps=self.locs)

    def addRightFrm(self):
        rightFrm = ttk.Frame(self)
        rightFrm.grid(row=0, column=3, rowspan=2, sticky="nsew")
        rightFrm.columnconfigure(0, weight=1)
        rightFrm.rowconfigure(0, weight=1)

        validationNN = self.register(validateNN)
        validationPI = self.register(validatePI)
        validationCE = self.register(validateCE)

        i = 1
        envFrm = LocLabelFrame(rightFrm, locKey="envFrmLabel", locFunc=self.getLocStr, allLLF=self.locs)
        envFrm.grid(row=i, column=0, sticky="nsew")
        envFrm.columnconfigure(0, weight=1)
        i += 1

        j = 0
        self.inAtmos = LocLabelCheck(
            parent=envFrm,
            labelLocKey="atmosLabel",
            row=j,
            col=0,
            locFunc=self.getLocStr,
            allLC=self.locs,
            columnspan=3,
        )
        self.inAtmos.trace_add("write", self.ambCallback)

        j += 1
        self.ambP = Loc3Input(
            parent=envFrm,
            row=j,
            col=0,
            labelLocKey="ambPresLabel",
            unitText="kPa",
            default="101.325",
            validation=validationNN,
            locFunc=self.getLocStr,
            allInputs=self.locs,
        )
        j += 1
        self.ambRho = Loc3Input(
            parent=envFrm,
            row=j,
            col=0,
            labelLocKey="ambRhoLabel",
            unitText="kg/m³",
            default="1.204",
            validation=validationNN,
            locFunc=self.getLocStr,
            allInputs=self.locs,
        )

        j += 1
        self.ambGam = Loc3Input(
            parent=envFrm,
            row=j,
            col=0,
            labelLocKey="ambGamLabel",
            default="1.400",
            validation=validationNN,
            locFunc=self.getLocStr,
            allInputs=self.locs,
        )

        solFrm = LocLabelFrame(
            rightFrm,
            locKey="solFrmLabel",
            locFunc=self.getLocStr,
            allLLF=self.locs,
        )
        solFrm.grid(row=i, column=0, sticky="nsew")
        i += 1
        solFrm.columnconfigure(0, weight=1)

        # noinspection SpellCheckingInspection
        self.dropSoln = LocDropdown(
            parent=solFrm,
            strObjDict={soln: soln for soln in (SOL_LAGRANGE, SOL_PIDDUCK, SOL_MAMONTOV)},
            locFunc=self.getLocStr,
            dropdowns=self.locs,
            descLabelKey="solFrmLabel",
        )
        self.dropSoln.grid(row=0, column=0, columnspan=2, sticky="nsew", padx=2, pady=2)

        opFrm = LocLabelFrame(
            rightFrm,
            locKey="opFrmLabel",
            locFunc=self.getLocStr,
            allLLF=self.locs,
        )
        opFrm.grid(row=4, column=0, sticky="nsew")
        opFrm.columnconfigure(0, weight=1)

        consFrm = LocLabelFrame(
            opFrm,
            locKey="consFrmLabel",
            style="SubLabelFrame.TLabelframe",
            locFunc=self.getLocStr,
            allLLF=self.locs,
        )
        consFrm.grid(row=i, column=0, columnspan=3, sticky="nsew", padx=2, pady=2)
        consFrm.columnconfigure(0, weight=1)
        i += 1

        j = 0
        self.solve_W_Lg = LocLabelCheck(
            parent=consFrm,
            row=j,
            col=0,
            columnspan=3,
            default=0,
            labelLocKey="consButton",
            tooltipLocKey="useConsText",
            locFunc=self.getLocStr,
            allLC=self.locs,
        )

        self.solve_W_Lg.trace_add("write", self.ctrlCallback)

        j += 1

        self.lock_Lg = LocLabelCheck(
            parent=consFrm,
            row=j,
            col=0,
            columnspan=3,
            default=0,
            labelLocKey="lockButton",
            tooltipLocKey="lockText",
            locFunc=self.getLocStr,
            allLC=self.locs,
        )
        self.lock_Lg.trace_add("write", self.ctrlCallback)
        j += 1

        self.opt_lf = LocLabelCheck(
            parent=consFrm,
            row=j,
            col=0,
            columnspan=3,
            default=0,
            labelLocKey="minTVButton",
            tooltipLocKey="optLFText",
            locFunc=self.getLocStr,
            allLC=self.locs,
        )
        self.opt_lf.trace_add("write", self.ctrlCallback)
        j += 1

        self.vTgt = Loc3Input(
            parent=consFrm,
            row=j,
            col=0,
            labelLocKey="vTgtLabel",
            unitText="m/s",
            default="1000.0",
            validation=validationNN,
            locFunc=self.getLocStr,
            allInputs=self.locs,
        )
        j += 1
        self.pTgt = Loc3Input(
            parent=consFrm,
            row=j,
            col=0,
            labelLocKey="pTgtLabel",
            unitText="MPa",
            default="350.0",
            validation=validationNN,
            tooltipLocKey="pTgtText",
            locFunc=self.getLocStr,
            allInputs=self.locs,
        )
        j += 1

        self.pControl = LocDropdown(
            parent=consFrm,
            strObjDict={" ": " "},
            locFunc=self.getLocStr,
            dropdowns=self.locs,
            descLabelKey="Pressure Constraint",
        )
        self.pControl.grid(row=j, column=0, columnspan=3, sticky="nsew", padx=2, pady=2)
        j += 1

        self.minWeb = Loc3Input(
            parent=consFrm,
            row=j,
            col=0,
            labelLocKey="minWebLabel",
            unitText="μm",
            default="1.0",
            validation=validationNN,
            color="red",
            locFunc=self.getLocStr,
            allInputs=self.locs,
        )
        j += 1
        # noinspection SpellCheckingInspection
        self.lgmax = Loc3Input(
            parent=consFrm,
            row=j,
            col=0,
            labelLocKey="maxLgLabel",
            unitText="m",
            default="100.0",
            validation=validationNN,
            color="red",
            locFunc=self.getLocStr,
            allInputs=self.locs,
        )

        i += 1
        sampleFrm = LocLabelFrame(
            opFrm,
            locKey="sampleFrmLabel",
            style="SubLabelFrame.TLabelframe",
            tooltipLocKey="sampText",
            locFunc=self.getLocStr,
            allLLF=self.locs,
        )
        sampleFrm.grid(row=i, column=0, columnspan=2, sticky="nsew", padx=2, pady=2)
        sampleFrm.columnconfigure(0, weight=1)
        self.dropDomain = LocDropdown(
            parent=sampleFrm,
            strObjDict={domain: domain for domain in (DOMAIN_TIME, DOMAIN_LEN)},
            locFunc=self.getLocStr,
            dropdowns=self.locs,
            descLabelKey="sampleFrmLabel",
        )
        j = 0
        self.dropDomain.grid(row=j, column=0, columnspan=2, sticky="nsew", padx=2, pady=2)
        j += 1
        self.step = Loc2Input(
            parent=sampleFrm,
            row=j,
            col=0,
            labelLocKey="stepLabel",
            default="33",
            validation=validationNN,
            formatter=formatIntInput,
            locFunc=self.getLocStr,
            allInputs=self.locs,
        )
        i += 1

        self.accExp = Loc2Input(
            parent=opFrm,
            row=i,
            col=0,
            labelLocKey="-log10(ε)",
            default="3",
            validation=validationPI,
            formatter=formatIntInput,
            color="red",
            tooltipLocKey="tolText",
            locFunc=self.getLocStr,
            allInputs=self.locs,
        )
        i += 1

        self.calcButton = ttk.Button(opFrm, text=self.getLocStr("calcLabel"), command=self.onCalculate)
        self.calcButton.grid(row=i, column=0, columnspan=3, sticky="nsew", padx=2, pady=2)
        # self.root.bind("<Return>", lambda *_: self.onCalculate())
        opFrm.rowconfigure(i, weight=1)
        i += 1
        self.progress = IntVar()
        ttk.Progressbar(opFrm, maximum=100, variable=self.progress).grid(
            row=i, column=0, columnspan=3, sticky="nsew", padx=2, pady=2
        )
        self.calcButtonTip = StringVar(value=self.getLocStr("calcButtonText"))
        CreateToolTip(self.calcButton, self.calcButtonTip)
        i += 1

        ctrlFrm = LocLabelFrame(rightFrm, locKey="guideCtrlFrmLabel", locFunc=self.getLocStr, allLLF=self.locs)
        ctrlFrm.grid(row=i, column=0, columnspan=3, sticky="nsew")
        ctrlFrm.columnconfigure(0, weight=1)
        i += 1

        self.guideMinLF, self.guideMaxLF, self.guideStepLF, self.guideMinCMR, self.guideMaxCMR, self.guideStepCMR = (
            Loc3Input(
                ctrlFrm,
                locFunc=self.getLocStr,
                allInputs=self.locs,
                labelLocKey=locKey,
                default=default,
                unitText=unit,
                validation=validation,
                row=j,
                col=0,
            )
            for j, (locKey, default, unit, validation) in enumerate(
                (
                    ("minLFLabel", "10.0", "%", validationCE),
                    ("maxLFLabel", "90.0", "%", validationCE),
                    ("stepLFLabel", "5.0", "%", validationCE),
                    ("minCMRLabel", "0.25", "", validationNN),
                    ("maxCMRLabel", "1.75", "", validationNN),
                    ("stepCMRLabel", "0.05", "", validationNN),
                )
            )
        )
        self.guideButton = ttk.Button(ctrlFrm, text=self.getLocStr("guideLabel"), command=self.onGuide)
        # self.root.bind("<space>", lambda *_: self.onGuide())
        self.guideButton.grid(row=6, column=0, columnspan=3, sticky="nsew", padx=2, pady=2)

    def generateKwargs(self):
        constrain = self.solve_W_Lg.get() == 1
        lock = self.lock_Lg.get() == 1
        optimize = self.opt_lf.get() == 1
        debug = self.debug.get() == 1
        atmosphere = self.inAtmos.get() == 1
        autofrettage = self.isAf.get() == 1

        gunType = self.typeOptn.getObj()

        if self.prop is None:
            raise ValueError("Invalid propellant.")

        if self.useCv.getObj() == USE_CV:
            chamberVolume = float(self.cvL.get()) * 1e-3
        else:
            chamberVolume = float(self.chgkg.get()) / self.prop.rho_p / float(self.ldf.get()) * 100

        chambrage = float(self.clr.get())
        chargeMass = float(self.chgkg.get())
        caliber = float(self.calmm.get()) * 1e-3

        gunLength = float(self.tblmm.get()) * 1e-3
        loadFraction = float(self.ldf.get()) * 1e-2

        self.kwargs = {
            "opt": optimize,
            "con": constrain,
            "deb": debug,
            "lock": lock,
            "typ": gunType,
            "dom": self.dropDomain.getObj(),
            "sol": self.dropSoln.getObj(),
            "control": self.pControl.getObj(),
            "structuralMaterial": self.dropMat.getObj().createMaterialAtTemp(self.dropMatTemp.getObj()),
            "structuralSafetyFactor": float(self.ssf.get()),
            "caliber": caliber,
            "shotMass": float(self.shtkg.get()),
            "propellant": self.prop,
            "grainSize": float(self.arcmm.get()) * 1e-3,
            "chargeMass": chargeMass,
            "chargeMassRatio": float(self.chgkg.get()) / float(self.shtkg.get()),
            "chamberVolume": chamberVolume,
            "startPressure": float(self.stpMPa.get()) * 1e6,
            "lengthGun": gunLength,
            "chambrage": chambrage,  # chamber expansion
            "nozzleExpansion": float(self.nozzExp.get()),  # nozzle expansion
            "nozzleEfficiency": float(self.nozzEff.get()) * 1e-2,  # nozzle efficiency
            "dragCoefficient": float(self.dgc.get()) * 1e-2,  # drag coefficient
            "designPressure": float(self.pTgt.get()) * 1e6,  # design pressure
            "designVelocity": float(self.vTgt.get()),  # design velocity
            "tol": 10 ** -int(self.accExp.get()),
            "minWeb": 1e-6 * float(self.minWeb.get()),
            "maxLength": float(self.lgmax.get()),
            "loadFraction": loadFraction,
            "step": int(self.step.get()),
            "autofrettage": autofrettage,
            "knownBore": lock,
            "minCMR": float(self.guideMinCMR.get()),
            "maxCMR": float(self.guideMaxCMR.get()),
            "stepCMR": float(self.guideStepCMR.get()),
            "minLF": float(self.guideMinLF.get()) * 1e-2,
            "maxLF": float(self.guideMaxLF.get()) * 1e-2,
            "stepLF": float(self.guideStepLF.get()) * 1e-2,
        }

        if atmosphere:
            self.kwargs.update(
                {
                    "ambientP": float(self.ambP.get()) * 1e3,
                    "ambientRho": float(self.ambRho.get()),
                    "ambientGamma": float(self.ambGam.get()),
                }
            )
        else:
            self.kwargs.update({"ambientP": 0, "ambientRho": 0, "ambientGamma": 1})

    def onGuide(self):
        if self.process or self.guideProcess:
            return

        self.focus()  # remove focus to force widget entry validation
        self.update()  # and wait for the event to update.

        self.guideProcess = None
        try:
            self.generateKwargs()
            self.guideProcess = Process(
                target=guide, args=(self.guideJobQueue, self.progressQueue, self.logQueue, self.kwargs)
            )
            self.guideProcess.start()

        except Exception:
            self.guide = None
            exc_type, exc_value, exc_traceback = sys.exc_info()
            rootLogger.error("exception when dispatching calculation:")
            rootLogger.error("".join(traceback.format_exception(exc_type, exc_value, exc_traceback)))

            self.updateGuideGraph()
        else:
            for loc in self.locs:
                try:
                    loc.inhibit()
                except AttributeError:
                    pass

            self.calcButton.config(state="disabled")
            self.guideButton.config(state="disabled")

    def onCalculate(self):
        if self.process or self.guideProcess:
            return
        self.focus()  # remove focus to force widget entry validation
        self.update()  # and wait for the event to update.

        self.process = None
        try:
            self.generateKwargs()
            self.process = Process(
                target=calculate, args=(self.jobQueue, self.progressQueue, self.logQueue, self.kwargs)
            )

            self.process.start()

        except Exception:
            exc_type, exc_value, exc_traceback = sys.exc_info()
            rootLogger.error("exception when dispatching calculation:")
            rootLogger.error("".join(traceback.format_exception(exc_type, exc_value, exc_traceback)))

            self.gun = None
            self.gunResult = None

            self.updateTable()
            self.updateFigPlot()
            self.updateAuxPlot()
        else:
            for loc in self.locs:
                try:
                    loc.inhibit()
                except AttributeError:
                    pass

            self.calcButton.config(state="disabled")
            self.guideButton.config(state="disabled")

    def getValue(self):
        try:
            self.kwargs, self.gun, self.gunResult = self.jobQueue.get_nowait()
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
        trueLF = kwargs["loadFraction"]

        boreS = pi * 0.25 * caliber**2
        breechS = boreS * chambrage

        sigfig = int(-log10(kwargs["tol"])) + 1
        gun = self.gun

        if gun:
            if constrain:
                # noinspection SpellCheckingInspection
                webmm = roundSig(kwargs["grainSize"] * 1e3, n=sigfig)
                self.arcmm.set(webmm)

                if not lock:
                    # noinspection SpellCheckingInspection
                    lgmm = roundSig(kwargs["lengthGun"] * 1e3, n=sigfig)
                    self.tblmm.set(lgmm)

                if optimize:
                    if self.useCv.getObj() == USE_CV:
                        self.cvL.set(roundSig(kwargs["chamberVolume"] * 1e3, n=sigfig))
                    else:
                        self.ldf.set(roundSig(kwargs["loadFraction"] * 100, n=sigfig))  # corrected "bulk" load fraction

            self.ld.set(toSI(trueLF * self.prop.rho_p, useSN=False, unit="kg/m³"))
            # true load density

            ps = self.gunResult.readTableData(POINT_PEAK_SHOT).shotPressure

            eta_t, eta_b, eta_p = self.gunResult.getEff()

            self.te.set(str(round(eta_t * 100, 2)) + " %")
            self.be.set(str(round(eta_b * 100, 2)) + " %")
            self.pe.set(str(round(eta_p * 100, 2)) + " %")

            self.va.set(toSI(gun.v_j, unit="m/s"))

            cartridge_len = kwargs["chamberVolume"] / breechS  # is equivalent to chamber length

            self.lx.set(
                (
                    toSI(kwargs["lengthGun"] / caliber, unit="Cal"),
                    toSI((kwargs["lengthGun"] + cartridge_len / kwargs["chambrage"]) / kwargs["caliber"], unit="Cal"),
                )
            )

            self.ammo.set(toSI(cartridge_len, unit="m"))
            self.pa.set(toSI(ps * gun.S / gun.m, unit="m/s²"))
            self.name.set(self.getDescriptive())

            tube_mass = self.gunResult.tubeMass
            self.gm.set(formatMass(tube_mass))

            peakAverageEntry = self.gunResult.readTableData(POINT_PEAK_AVG)
            peakBreechEntry = self.gunResult.readTableData(POINT_PEAK_BREECH)

            self.pp.set(
                (toSI(peakAverageEntry.avgPressure, unit="Pa"), toSI(peakBreechEntry.breechPressure, unit="Pa"))
            )

            muzzleEntry = self.gunResult.readTableData(POINT_EXIT)

            self.mv.set(toSI(muzzleEntry.velocity, unit="m/s"))

            # p = self.tabParent.index(self.tabParent.select())  # find the currently active tab
            # if p == 2 or p == 3:  # if the guidance/logging pane is currently selected.
            #     self.tabParent.select(0)  # goto the graphing pane

        else:  # calculation results in some error
            # self.tabParent.select(3)  # goto the error pane
            pass

        self.updateTable()
        self.updateFigPlot()
        self.updateAuxPlot()
        self.calcButton.config(state="normal")
        self.guideButton.config(state="normal")

        for loc in self.locs:
            try:
                loc.disinhibit()
            except AttributeError:
                pass

    def getGuide(self):
        try:
            self.guide = self.guideJobQueue.get_nowait()
        except Empty:
            return

        self.guideProcess = None
        self.updateGuideGraph()
        self.calcButton.config(state="normal")
        self.guideButton.config(state="normal")

        for loc in self.locs:
            try:
                loc.disinhibit()
            except AttributeError:
                pass

    def addSpecFrm(self):
        specFrm = LocLabelFrame(self, locKey="specFrmLabel", locFunc=self.getLocStr, allLLF=self.locs)
        specFrm.grid(row=0, column=2, rowspan=2, sticky="nsew")
        specFrm.columnconfigure(0, weight=1)

        # validation
        validationNN = self.register(validateNN)
        validationCE = self.register(validateCE)

        i = 0
        # noinspection SpellCheckingInspection
        self.typeOptn = LocDropdown(
            parent=specFrm,
            strObjDict={gun_type: gun_type for gun_type in (CONVENTIONAL, RECOILLESS)},
            locFunc=self.getLocStr,
            dropdowns=self.locs,
            descLabelKey="typeLabel",
        )

        self.typeOptn.grid(row=i, column=0, sticky="nsew", padx=2, pady=2, columnspan=3)
        i += 1
        # noinspection SpellCheckingInspection
        self.calmm = Loc3Input(
            parent=specFrm,
            row=i,
            labelLocKey="calLabel",
            unitText="mm",
            default="50.0",
            validation=validationNN,
            locFunc=self.getLocStr,
            allInputs=self.locs,
        )
        i += 1
        # noinspection SpellCheckingInspection
        self.tblmm = Loc3Input(
            parent=specFrm,
            row=i,
            labelLocKey="tblLabel",
            unitText="mm",
            default="3500.0",
            validation=validationNN,
            locFunc=self.getLocStr,
            allInputs=self.locs,
        )
        i += 1
        # noinspection SpellCheckingInspection
        self.shtkg = Loc3Input(
            parent=specFrm,
            row=i,
            labelLocKey="shtLabel",
            unitText="kg",
            default="1.0",
            validation=validationNN,
            locFunc=self.getLocStr,
            allInputs=self.locs,
        )
        i += 1
        # noinspection SpellCheckingInspection
        self.chgkg = Loc3Input(
            parent=specFrm,
            row=i,
            labelLocKey="chgLabel",
            unitText="kg",
            default="0.5",
            validation=validationNN,
            tooltipLocKey="chgText",
            locFunc=self.getLocStr,
            allInputs=self.locs,
        )
        i += 1

        self.propTabParent = ttk.Notebook(specFrm, padding=0)
        self.propTabParent.grid(row=i, column=0, columnspan=3, sticky="nsew")
        specFrm.rowconfigure(i, weight=1)

        self.propTabParent.columnconfigure(0, weight=1)
        self.propTabParent.rowconfigure(0, weight=1)

        propFrm = LocLabelFrame(
            self.propTabParent,
            locKey="propFrmLabel",
            style="SubLabelFrame.TLabelframe",
            locFunc=self.getLocStr,
            tooltipLocKey="specsText",
            allLLF=self.locs,
        )
        self.propFrm = propFrm
        propFrm.grid(row=0, column=0, sticky="nsew", padx=2, pady=2)
        propFrm.rowconfigure(1, weight=1)
        propFrm.columnconfigure(0, weight=1)

        self.dropProp = LocDropdown(
            parent=propFrm,
            strObjDict=GrainComp.readFile(
                resolvepath("ballistics/resource/propellants.csv")
            ),  # dict of composition.name (string) -> composition (object),
            locFunc=self.getLocStr,
            dropdowns=self.locs,
            descLabelKey="propFrmLabel",
        )
        self.dropProp.grid(row=0, column=0, columnspan=2, sticky="nsew", padx=2, pady=2)

        specScroll = ttk.Scrollbar(propFrm, orient="vertical")
        specScroll.grid(row=1, rowspan=2, column=1, sticky="nsew")
        specHScroll = ttk.Scrollbar(propFrm, orient="horizontal")
        specHScroll.grid(row=2, column=0, sticky="nsew")

        self.specs = Text(
            propFrm,
            wrap="word",
            height=10,
            width=30,
            yscrollcommand=specScroll.set,
            xscrollcommand=specHScroll.set,
            font=(FONTNAME, FONTSIZE),
        )
        self.specs.grid(row=1, column=0, sticky="nsew")
        specScroll.config(command=self.specs.yview)
        specHScroll.config(command=self.specs.xview)

        self.forceUpdOnThemeWidget.append(self.specs)

        grainFrm = LocLabelFrame(
            self.propTabParent,
            locKey="grainFrmLabel",
            style="SubLabelFrame.TLabelframe",
            locFunc=self.getLocStr,
            allLLF=self.locs,
        )
        self.grainFrm = grainFrm
        grainFrm.grid(row=i, column=0, columnspan=3, sticky="nsew", padx=2, pady=2)
        grainFrm.columnconfigure(0, weight=1)
        grainFrm.rowconfigure(0, weight=1)
        self.propTabParent.add(propFrm, text=self.getLocStr("propFrmLabel"))
        self.propTabParent.add(grainFrm, text=self.getLocStr("grainFrmLabel"))

        self.propTabParent.enable_traversal()

        geomPlotFrm = LocLabelFrame(
            grainFrm,
            locKey="σ(Z)",
            style="SubLabelFrame.TLabelframe",
            tooltipLocKey="geomPlotText",
            locFunc=self.getLocStr,
            allLLF=self.locs,
        )

        j = 0
        geomPlotFrm.grid(row=j, column=0, columnspan=3, sticky="nsew", padx=2, pady=2)

        self.geomPlotFrm = geomPlotFrm
        self.dropGeom = LocDropdown(
            parent=grainFrm,
            strObjDict=GEOMETRIES,
            locFunc=self.getLocStr,
            dropdowns=self.locs,
            descLabelKey="grainFrmLabel",
        )

        j += 1
        self.dropGeom.grid(row=j, column=0, columnspan=3, sticky="nsew", padx=2, pady=2)
        j += 1
        # noinspection SpellCheckingInspection
        self.arcmm = Loc3Input(
            parent=grainFrm,
            row=j,
            labelLocKey="",
            descLabelKey="Web",
            unitText="mm",
            default="1.0",
            validation=validationNN,
            tooltipLocKey="",
            locFunc=self.getLocStr,
            allInputs=self.locs,
        )
        j += 1

        self.grdR = Loc3Input(
            parent=grainFrm,
            row=j,
            labelLocKey="",
            descLabelKey="1/α",
            unitText="x",
            default="1.0",
            validation=validationNN,
            tooltipLocKey="",
            locFunc=self.getLocStr,
            allInputs=self.locs,
        )
        j += 1
        self.grlR = Loc3Input(
            parent=grainFrm,
            row=j,
            labelLocKey="",
            descLabelKey="1/β",
            unitText="x",
            default="2.5",
            validation=validationNN,
            tooltipLocKey="",
            locFunc=self.getLocStr,
            allInputs=self.locs,
        )

        i += 1
        # noinspection SpellCheckingInspection
        self.useCv = LocDropdown(
            parent=specFrm,
            strObjDict={USE_LF: USE_LF, USE_CV: USE_CV},
            locFunc=self.getLocStr,
            dropdowns=self.locs,
            descLabelKey="cvlfLabel",
        )
        self.useCv.grid(row=i, column=0, columnspan=3, sticky="nsew", padx=2, pady=2)

        i += 1
        self.ldf = Loc3Input(
            parent=specFrm,
            row=i,
            labelLocKey="ldfLabel",
            unitText="%",
            default="50.0",
            validation=validationCE,
            tooltipLocKey="ldfText",
            locFunc=self.getLocStr,
            allInputs=self.locs,
        )
        i += 1
        self.cvL = Loc3Input(
            parent=specFrm,
            row=i,
            labelLocKey="cvLabel",
            unitText="L",
            default="20.0",
            validation=validationNN,
            locFunc=self.getLocStr,
            allInputs=self.locs,
        )

        self.useCv.trace_add("write", self.cvlfCallback)
        for widget in (self.cvL, self.ldf, self.chgkg, self.accExp):
            widget.trace_add("write", self.cvlfConsistencyCallback)

        i += 1
        self.clr = Loc3Input(
            parent=specFrm,
            row=i,
            labelLocKey="clrLabel",
            unitText="x",
            default="1.5",
            validation=validationNN,
            tooltipLocKey="clrText",
            locFunc=self.getLocStr,
            allInputs=self.locs,
        )
        i += 1
        self.dgc = Loc3Input(
            parent=specFrm,
            row=i,
            labelLocKey="dgcLabel",
            unitText="%",
            default="3.0",
            validation=validationCE,
            tooltipLocKey="dgcText",
            locFunc=self.getLocStr,
            allInputs=self.locs,
        )
        i += 1
        self.stpMPa = Loc3Input(
            parent=specFrm,
            row=i,
            labelLocKey="stpLabel",
            unitText="MPa",
            default="30.0",
            validation=validationNN,
            tooltipLocKey="stpText",
            locFunc=self.getLocStr,
            allInputs=self.locs,
        )
        i += 1

        # noinspection SpellCheckingInspection
        self.nozzExp = Loc3Input(
            parent=specFrm,
            row=i,
            labelLocKey="nozzExpLabel",
            unitText="x",
            default="4.0",
            validation=validationNN,
            tooltipLocKey="nozzExpText",
            locFunc=self.getLocStr,
            allInputs=self.locs,
        )
        i += 1
        self.nozzEff = Loc3Input(
            parent=specFrm,
            row=i,
            labelLocKey="nozzEffLabel",
            unitText="%",
            default="92.0",
            validation=validationCE,
            tooltipLocKey="nozzEffText",
            locFunc=self.getLocStr,
            allInputs=self.locs,
        )
        i += 1

        mecFrm = LocLabelFrame(specFrm, locKey="matFrmLabel", locFunc=self.getLocStr, allLLF=self.locs)
        mecFrm.grid(row=i, column=0, sticky="nsew", columnspan=3)
        mecFrm.columnconfigure(0, weight=1)
        i += 1

        j = 0
        self.dropMat = LocDropdown(
            parent=mecFrm, strObjDict=MATERIALS, locFunc=self.getLocStr, dropdowns=self.locs, descLabelKey="matFrmLabel"
        )
        self.dropMat.grid(row=j, column=0, columnspan=2, sticky="nsew", padx=2, pady=2)
        self.dropMat.trace_add("write", lambda *args: self.dropMatTemp.reset(self.dropMat.getObj().getTdict()))
        j += 1
        self.dropMatTemp = LocDropdown(
            parent=mecFrm,
            strObjDict=self.dropMat.getObj().getTdict(),
            locFunc=self.getLocStr,
            dropdowns=self.locs,
            descLabelKey="matFrmTempLabel",
        )
        self.dropMatTemp.grid(row=j, column=0, columnspan=2, sticky="nsew", padx=2, pady=2)
        j += 1

        self.ssf = Loc2Input(
            parent=mecFrm,
            row=j,
            col=0,
            labelLocKey="sffLabel",
            default="1.35",
            validation=validationNN,
            locFunc=self.getLocStr,
            allInputs=self.locs,
        )
        j += 1
        self.isAf = LocLabelCheck(
            parent=mecFrm,
            labelLocKey="afLabel",
            row=j,
            col=0,
            locFunc=self.getLocStr,
            allLC=self.locs,
            columnspan=2,
        )
        j += 1

        i += 1
        specFrm.rowconfigure(i, weight=1)

        self.dropProp.trace_add("write", self.updateSpec)
        self.dropGeom.trace_add("write", self.updateGeom)

        self.grdR.trace_add("write", self.callback)
        self.grlR.trace_add("write", self.callback)
        self.arcmm.trace_add("write", self.callback)

        self.typeOptn.trace_add("write", self.typeCallback)

    def addGeomPlot(self):
        with mpl.rc_context(CONTEXT):
            fig = Figure(dpi=96, layout="constrained")
            self.geomFig = fig
            self.geomAx = fig.add_subplot(111)

            self.geomCanvas = FigureCanvasTkAgg(fig, master=self.geomPlotFrm)
            self.geomCanvas.draw_idle()
            self.geomCanvas.get_tk_widget().place(relheight=1, relwidth=1)

    def addGuidePlot(self):
        plotFrm = LocLabelFrame(self.guideTab, locKey="guideFrmLabel", locFunc=self.getLocStr, allLLF=self.locs)
        plotFrm.grid(row=0, column=0, sticky="nsew")

        with mpl.rc_context(CONTEXT):
            fig = Figure(dpi=96, layout="constrained")
            self.guideFig = fig
            self.guideAx = fig.add_subplot(111)
            self.guideCanvas = FigureCanvasTkAgg(fig, master=plotFrm)
            self.guideCanvas.draw_idle()
            self.guideCanvas.get_tk_widget().place(relheight=1, relwidth=1)

        controlFrm = Frame(self.guideTab)
        controlFrm.grid(row=1, column=0, sticky="nsew")

        self.guideIndex = IntVar(value=3)
        self.guideIndex.trace_add("write", callback=self.guideCallback)

        for i in range(3):
            controlFrm.columnconfigure(i, weight=1)
        self.guidePlotTravel = ttk.Radiobutton(
            controlFrm, text=self.getLocStr("guidePlotTravel"), variable=self.guideIndex, value=3
        )
        self.guidePlotTravel.grid(row=0, column=0, sticky="nsew")

        self.guidePlotVolume = ttk.Radiobutton(
            controlFrm, text=self.getLocStr("guidePlotVolume"), variable=self.guideIndex, value=4
        )
        self.guidePlotVolume.grid(row=0, column=1, sticky="nsew")

        self.guidePlotBurnout = ttk.Radiobutton(
            controlFrm, text=self.getLocStr("guidePlotBurnout"), variable=self.guideIndex, value=5
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

        self.tabParent.add(self.plotTab, text=self.getLocStr("plotTab"))
        self.tabParent.add(self.tableTab, text=self.getLocStr("tableTab"))
        self.tabParent.add(self.guideTab, text=self.getLocStr("guideTab"))
        self.tabParent.add(self.errorTab, text=self.getLocStr("errorTab"))

        self.tabParent.enable_traversal()

        self.addTblFrm()
        self.addErrFrm()

        self.addFigPlot()

    def addErrFrm(self):
        errorFrm = LocLabelFrame(self.errorTab, locKey="errFrmLabel", locFunc=self.getLocStr, allLLF=self.locs)
        errorFrm.grid(row=0, column=0, columnspan=3, sticky="nsew")
        errorFrm.columnconfigure(0, weight=1)
        errorFrm.rowconfigure(0, weight=1)

        errScroll = ttk.Scrollbar(errorFrm, orient="vertical")
        errScroll.grid(row=0, column=1, sticky="nsew")
        self.errorText = Text(
            errorFrm,
            yscrollcommand=errScroll.set,
            wrap="word",
            height=6,
            width=0,
            font=(FONTNAME, FONTSIZE),
            state="disabled",
        )
        self.errorText.grid(row=0, column=0, sticky="nsew")

        self.errorText.tag_configure(str(logging.DEBUG), foreground="tan")
        self.errorText.tag_configure(str(logging.WARNING), foreground="orange")
        self.errorText.tag_configure(str(logging.ERROR), foreground="orangered")
        self.errorText.tag_configure(str(logging.CRITICAL), foreground="red")

        errScroll.config(command=self.errorText.yview)
        self.forceUpdOnThemeWidget.append(self.errorText)

    def addFigPlot(self):
        plotFrm = LocLabelFrame(
            self.plotTab,
            locKey="plotFrmLabel",
            locFunc=self.getLocStr,
            tooltipLocKey="plotText",
            allLLF=self.locs,
        )
        plotFrm.grid(row=0, column=0, sticky="nsew")

        for i in range(5):
            plotFrm.columnconfigure(i, weight=1)

        checks = []

        j = 1
        k = 0
        self.plotAvgP = LocLabelCheck(
            parent=plotFrm, labelLocKey="plotAvgP", row=j, col=k, locFunc=self.getLocStr, allLC=checks
        )
        k += 1
        self.plotBaseP = LocLabelCheck(
            parent=plotFrm, labelLocKey="plotBaseP", row=j, col=k, locFunc=self.getLocStr, allLC=checks
        )
        k += 1
        self.plotBreechP = LocLabelCheck(
            parent=plotFrm, labelLocKey="plotBreechP", row=j, col=k, locFunc=self.getLocStr, allLC=checks
        )
        k += 1
        self.plotStagP = LocLabelCheck(
            parent=plotFrm, labelLocKey="plotStagP", row=j, col=k, locFunc=self.getLocStr, allLC=checks
        )

        j += 1
        k = 0
        self.plotVel = LocLabelCheck(
            parent=plotFrm, labelLocKey="plotVel", row=j, col=k, locFunc=self.getLocStr, allLC=checks
        )
        k += 1
        self.plotNozzleV = LocLabelCheck(
            parent=plotFrm, labelLocKey="plotNozzleV", row=j, col=k, locFunc=self.getLocStr, allLC=checks
        )
        k += 1
        self.plotBurnup = LocLabelCheck(
            parent=plotFrm, labelLocKey="plotBurnup", row=j, col=k, locFunc=self.getLocStr, allLC=checks
        )
        k += 1
        self.plotEta = LocLabelCheck(
            parent=plotFrm, labelLocKey="plotEta", row=j, col=k, locFunc=self.getLocStr, allLC=checks
        )

        for check in checks:
            check.trace_add("write", self.updateFigPlot)

        plotFrm.columnconfigure(0, weight=1)
        plotFrm.rowconfigure(0, weight=1)

        plotPlaceFrm = Frame(plotFrm)
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

        auxFrm = LocLabelFrame(
            self.plotTab, locKey="auxFrmLabel", locFunc=self.getLocStr, tooltipLocKey="auxText", allLLF=self.locs
        )
        auxFrm.grid(row=1, column=0, sticky="nsew")

        for i in range(2):
            auxFrm.columnconfigure(i, weight=1)

        auxChecks = []
        j = 1
        k = 0
        self.traceHull = LocLabelCheck(
            parent=auxFrm, row=j, col=k, labelLocKey="traceHull", locFunc=self.getLocStr, default=1, allLC=auxChecks
        )

        k += 1
        self.tracePress = LocLabelCheck(
            parent=auxFrm, row=j, col=k, labelLocKey="tracePress", locFunc=self.getLocStr, allLC=auxChecks
        )

        for check in auxChecks:
            check.trace_add("write", self.updateAuxPlot)

        auxFrm.columnconfigure(0, weight=1)
        auxFrm.rowconfigure(0, weight=1)

        auxPlaceFrm = Frame(auxFrm)
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
                    for tag, t, l, psi, v, Pb, P, Ps, T in self.gunResult.getRawTableData():
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
                    for tag, t, l, psi, v, vx, Px, P0, P, Ps, T, eta in self.gunResult.getRawTableData():
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

                if self.plotBreechP.get():
                    self.axP.plot(
                        xs,
                        Pbs,
                        c="xkcd:goldenrod",
                        label=(
                            self.getLocStr("figBreech")
                            if gunType == CONVENTIONAL
                            else (
                                self.getLocStr("figNozzleP") if gunType == RECOILLESS else self.getLocStr("figBleedP")
                            )
                        ),
                    )

                if gunType == RECOILLESS:
                    if self.plotStagP.get():
                        self.axP.plot(xs, P0s, "seagreen", label=self.getLocStr("figStagnation"))

                    if self.plotNozzleV.get():
                        self.axv.plot(xs, vxs, "royalblue", label=self.getLocStr("figNozzleV"))

                    if self.plotEta.get():
                        self.ax.plot(xs, etas, "crimson", label=self.getLocStr("figOutflow"))

                if self.plotAvgP.get():
                    self.axP.plot(xs, Pas, "tab:green", label=self.getLocStr("figAvgP"))

                if self.plotBaseP.get():
                    self.axP.plot(xs, Pss, "yellowgreen", label=self.getLocStr("figShotBase"))

                if gunType == CONVENTIONAL or gunType == RECOILLESS:
                    self.axP.axhline(
                        float(self.pTgt.get()), c="tab:green", linestyle=":", label=self.getLocStr("figTgtP")
                    )

                if self.plotVel.get():
                    self.axv.plot(xs, vs, "tab:blue", label=self.getLocStr("figShotVel"))
                self.axv.axhline(vTgt, c="tab:blue", linestyle=":", label=self.getLocStr("figTgtV"))

                if self.plotBurnup.get():
                    self.ax.plot(xs, psis, c="tab:red", label=self.getLocStr("figPsi"))

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
                    self.ax.set_xlabel(self.getLocStr("figTimeDomain"))
                elif dom == DOMAIN_LEN:
                    self.ax.set_xlabel(self.getLocStr("figLenDomain"))

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

            pTrace = self.gunResult.pressureTrace
            # noinspection SpellCheckingInspection
            cmap = mpl.colormaps["afmhot" + ("_r" if THEMES[self.themeNameVar.get()]  else "")]

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

            self.auxAx.set_xlabel(self.getLocStr("figAuxDomain"))

            self.auxAx.yaxis.label.set_color("tab:green")
            self.auxAx.set_ylabel("MPa")

            self.auxAxH.yaxis.set_ticks_position("right")
            self.auxAxH.yaxis.set_label_position("right")

            self.auxAxH.yaxis.label.set_color("tab:blue")
            self.auxAxH.set_ylabel("mm")

            HTrace = self.gunResult.outline

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
        tblFrm = LocLabelFrame(self.tableTab, locKey="tblFrmLabel", locFunc=self.getLocStr, allLLF=self.locs)
        tblFrm.grid(row=0, column=0, sticky="nsew")

        tblFrm.columnconfigure(0, weight=1)
        tblFrm.rowconfigure(0, weight=1)

        tblPlaceFrm = Frame(tblFrm)
        tblPlaceFrm.grid(row=0, column=0, sticky="nsew")
        # noinspection SpellCheckingInspection
        vertscroll = ttk.Scrollbar(tblFrm, orient="vertical")  # create a scrollbar
        vertscroll.grid(row=0, rowspan=2, column=1, sticky="nsew")
        # noinspection SpellCheckingInspection
        horzscroll = ttk.Scrollbar(tblFrm, orient="horizontal")
        horzscroll.grid(row=1, column=0, sticky="nsew")

        self.tv = ttk.Treeview(
            tblPlaceFrm, selectmode="browse", yscrollcommand=vertscroll.set, xscrollcommand=horzscroll.set
        )  # this set the nbr. of values
        self.tv.place(relwidth=1, relheight=1)

        vertscroll.configure(command=self.tv.yview)  # make it vertical
        horzscroll.configure(command=self.tv.xview)

    # noinspection PyUnusedLocal
    def updateSpec(self, *args):
        self.specs.config(state="normal")
        compo = self.dropProp.getObj()
        self.specs.delete("1.0", "end")
        t_Font = tkFont.Font(font=self.specs.cget("font"))
        width = self.specs.winfo_width() // t_Font.measure("m")
        # noinspection SpellCheckingInspection
        if compo.T_v:
            self.specs.insert(
                "end",
                "{:}: {:>4.0f} K {:}\n".format(self.getLocStr("TvDesc"), compo.T_v, self.getLocStr("isochorDesc")),
            )
        self.specs.insert("end", "{:}: {:>4.0f} kg/m³\n".format(self.getLocStr("densityDesc"), compo.rho_p))
        self.specs.insert("end", "{:}: {:>4.0f} kJ/kg\n".format(self.getLocStr("force"), compo.f / 1e3))
        isp = compo.getIsp()
        self.specs.insert(
            "end", "{:}: {:>4.0f} m/s {:>3.0f} s\n".format(self.getLocStr("vacISPDesc"), isp, isp / 9.805)
        )
        isp = compo.getIsp(50)
        self.specs.insert(
            "end",
            "{:}: {:>4.0f} m/s {:>3.0f} s\n{:}\n".format(
                self.getLocStr("atmISPDesc"), isp, isp / 9.805, self.getLocStr("pRatioDesc")
            ),
        )
        self.specs.insert("end", "{:}:\n".format(self.getLocStr("brDesc")))
        for p in (100e6, 200e6, 300e6):
            self.specs.insert(
                "end",
                "{:>12}".format(toSI(compo.getLBR(p), unit="m/s", dec=3))
                + " @ {:>12}\n".format(toSI(p, unit="Pa", dec=3)),
            )

        # self.specs.insert("end", "-" * width + "\n")
        for line in compo.desc.split(","):
            self.specs.insert("end", line + "\n")
        self.specs.config(state="disabled")
        # this updates the specification description

        self.callback()
        self.cvlfConsistencyCallback()  # update the chamber volume / load fraction with current data

    # noinspection PyUnusedLocal
    def updateGeom(self, *args):
        geom = self.dropGeom.getObj()
        if geom == SimpleGeometry.SPHERE:
            self.grlR.remove()
            self.grdR.remove()

        elif geom == SimpleGeometry.CYLINDER:
            self.grlR.restore()
            self.grdR.remove()

        else:
            self.grlR.restore()
            self.grdR.restore()

        if geom == SimpleGeometry.SPHERE:
            self.arcmm.reLocalize("diamLabel", "diaText")
            self.grlR.reLocalize("", "")
            self.grdR.reLocalize("", "")

        elif geom == SimpleGeometry.ROD:
            self.arcmm.reLocalize("widthLabel", "widthText")
            self.grlR.reLocalize("ltwLabel", "rodRText")
            self.grdR.reLocalize("htwLabel", "heightRText")

        elif geom == SimpleGeometry.CYLINDER:
            self.arcmm.reLocalize("diamLabel", "diaText")
            self.grlR.reLocalize("ltdLabel", "cylLRText")
            self.grdR.reLocalize("", "")

        else:
            self.arcmm.reLocalize("athLabel", "arcText")
            self.grlR.reLocalize("ltdLabel", "perfLRText")
            # noinspection SpellCheckingInspection
            self.grdR.reLocalize("pdtalLabel", "pDiaRText")

        self.callback()

    def updateGeomPlot(self):
        with mpl.rc_context(CONTEXT):
            N = 10
            prop = self.prop
            Zb = prop.Z_b
            self.geomAx.cla()
            if prop is not None:
                xs = [i / N for i in range(N)]
                ys = [prop.f_sigma_Z(x) for x in xs]

                if Zb > 1:
                    xs.extend((1, 1))
                    ys.extend(prop.f_ullim())

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
            tableData, errorData = (self.gunResult.getRawTableData(), self.gunResult.getRawErrorData())
        except AttributeError:
            gunType = self.typeOptn.getObj()
            tableData, errorData = [], []

        locTableData = []
        for i, line in enumerate(tableData):
            locTableData.append((self.getLocStr(line[0]), *line[1:]))

        if gunType == CONVENTIONAL:
            useSN = (False, False, False, True, False, False, False, False, True)
            units = (None, "s", "m", None, "m/s", "Pa", "Pa", "Pa", "K")
        elif gunType == RECOILLESS:
            useSN = (False, False, False, True, False, False, False, False, False, False, True, True)
            units = (None, "s", "m", None, "m/s", "m/s", "Pa", "Pa", "Pa", "Pa", "K", None)
        else:
            raise ValueError("unknown gun types")

        locTableData = dot_aligned(locTableData, units=units, useSN=useSN)
        errorData = dot_aligned(errorData, units=units, useSN=useSN)

        columnList = self.getLocStr("columnList")[gunType]
        self.tv["columns"] = columnList
        self.tv["show"] = "headings"

        self.tv.tag_configure(self.getLocStr(POINT_PEAK_STAG), foreground="#2e8b57")
        self.tv.tag_configure(self.getLocStr(POINT_PEAK_AVG), foreground="#2ca02c")
        self.tv.tag_configure(self.getLocStr(POINT_PEAK_BREECH), foreground="orange")
        self.tv.tag_configure(self.getLocStr(POINT_PEAK_SHOT), foreground="yellow green")
        self.tv.tag_configure(self.getLocStr(POINT_BURNOUT), foreground="red")
        self.tv.tag_configure(self.getLocStr(POINT_FRACTURE), foreground="brown")
        self.tv.tag_configure(self.getLocStr(POINT_EXIT), foreground="steel blue")
        self.tv.tag_configure(self.getLocStr(POINT_START), foreground="steel blue")
        self.tv.tag_configure(self.getLocStr(COMPUTE), foreground="tan")

        t_Font = tkFont.Font(family=FONTNAME, size=FONTSIZE)

        self.tv.tag_configure("monospace", font=t_Font)
        self.tv.tag_configure("error", font=t_Font, foreground="dim gray")

        # we use a fixed width font so any char will do
        fontWidth, _ = t_Font.measure("m"), t_Font.metrics("linespace")

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
            self.guideAx.cla()

            if self.guide and any(line for line in self.guide):
                loadDensities = list(list(value[0] for value in line) for line in self.guide)
                chargeMasses = list(list(value[1] for value in line) for line in self.guide)

                self.guideAx.set_xlabel(self.getLocStr("guideLDDomain"))
                self.guideAx.set_ylabel(self.getLocStr("guideCMDomain"))

                index = int(self.guideIndex.get())
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

                    self.guideAx.pcolormesh(
                        loadDensities,
                        chargeMasses,
                        values,
                        shading="nearest",
                        cmap="afmhot" + ("" if THEMES[self.themeNameVar.get()] else "_r"),
                        vmin=minVal,
                        vmax=maxVal,
                    )
                    threshold = 0.5 * (minVal + maxVal)
                    for line in self.guide:
                        for value in line:
                            loadDensity, chargeMass = value[0], value[1]
                            entry = value[index]
                            entry = entry * (1e3 if index == 4 else 1) if entry else None
                            self.guideAx.text(
                                loadDensity,
                                chargeMass,
                                (f"{entry:.3g}" if index != 5 else f"{entry:.3f}") if entry else "N/A",
                                color=bgc if (entry and entry < threshold) else fgc,
                                horizontalalignment="center",
                                verticalalignment="center",
                            )
                    if index == 3:
                        self.guideAx.set_title(self.getLocStr("guideTravelTitle"))
                    elif index == 4:
                        self.guideAx.set_title(self.getLocStr("guideBVTitle"))
                    else:
                        self.guideAx.set_title(self.getLocStr("guideBurnoutTitle"))

                else:
                    raise ValueError("unknown guidance diagram plotting control.")

                # self.guideCursor = Cursor(self.guideAx, useblit=True, color=fgc, linewidth=1)
            self.guideCanvas.draw_idle()

    # noinspection PyUnusedLocal
    def callback(self, *args):
        """
        updates the propellant object on write to the ratio entry fields
        and, on changing the propellant or geometrical specification.
        """
        geom = self.dropGeom.getObj()
        compo = self.dropProp.getObj()

        try:
            self.prop = Propellant(compo, geom, float(self.grdR.get()), float(self.grlR.get()))
            self.updateGeomPlot()

        except Exception:
            exc_type, exc_value, exc_traceback = sys.exc_info()
            self.prop = None
            rootLogger.error("exception in propellant callback:")
            rootLogger.error("".join(traceback.format_exception(exc_type, exc_value, exc_traceback)))

    # noinspection PyUnusedLocal
    def typeCallback(self, *args):
        gunType = self.typeOptn.getObj()

        if gunType == CONVENTIONAL:
            self.dropSoln.enable()
        else:
            self.dropSoln.setByObj(SOL_LAGRANGE)
            self.dropSoln.disable()

        if gunType == RECOILLESS:
            self.nozzExp.restore()
            self.nozzEff.restore()
            self.plotNozzleV.restore()

        else:
            self.nozzExp.remove()
            self.nozzEff.remove()
            self.plotNozzleV.remove()

        if gunType == CONVENTIONAL:
            self.plotBreechP.reLocalize("plotBreechP")
            self.plotStagP.remove()
            self.plotEta.remove()
            self.pControl.reset({point: point for point in (POINT_PEAK_AVG, POINT_PEAK_SHOT, POINT_PEAK_BREECH)})
        elif gunType == RECOILLESS:
            self.plotBreechP.reLocalize("plotNozzleP")
            self.plotStagP.restore()
            self.plotStagP.reLocalize("plotStagP")
            self.plotEta.restore()
            self.plotEta.reLocalize("plotEtaEsc")
            self.pControl.reset(
                {point: point for point in (POINT_PEAK_AVG, POINT_PEAK_SHOT, POINT_PEAK_STAG, POINT_PEAK_BREECH)}
            )
        else:
            raise ValueError("unknown gun type")

    # noinspection PyUnusedLocal
    def ctrlCallback(self, *args):
        if self.solve_W_Lg.get() == 0:
            self.vTgt.disable()
            self.pTgt.disable()

            self.opt_lf.disable()
            self.lock_Lg.disable()
            self.minWeb.disable()
            self.lgmax.disable()
            self.pControl.disable()
            self.tblmm.enable()

        else:
            if self.lock_Lg.get() == 1:
                self.vTgt.disable()
                self.opt_lf.disable()
                self.tblmm.enable()
            else:
                self.vTgt.enable()
                self.opt_lf.enable()
                self.tblmm.disable()

            if self.opt_lf.get() == 1:
                self.lock_Lg.disable()
            else:
                self.lock_Lg.enable()

            self.pTgt.enable()

            self.minWeb.enable()
            self.lgmax.enable()
            self.pControl.enable()

    # noinspection PyUnusedLocal, SpellCheckingInspection
    def cvlfConsistencyCallback(self, *args):
        try:
            sigfig = int(self.accExp.get()) + 1
            if self.useCv.getObj() == USE_CV:  # use Cv
                cv = float(self.cvL.get()) / 1e3
                w = float(self.chgkg.get())
                rho = self.prop.rho_p
                self.ldf.set(roundSig(w / cv / rho * 100, n=sigfig))

            else:  # using load fraction
                w = float(self.chgkg.get())
                lf = float(self.ldf.get()) / 100
                rho = self.prop.rho_p
                self.cvL.set(roundSig((w / rho / lf) * 1e3, n=sigfig))

        except (ZeroDivisionError, ValueError):
            return

    # noinspection PyUnusedLocal
    def ambCallback(self, *args):
        self.ambP.enable() if self.inAtmos.get() else self.ambP.disable()
        self.ambRho.enable() if self.inAtmos.get() else self.ambRho.disable()
        self.ambGam.enable() if self.inAtmos.get() else self.ambGam.disable()

    # noinspection PyUnusedLocal, SpellCheckingInspection
    def cvlfCallback(self, *args):
        useCv = self.useCv.getObj() == USE_CV

        self.ldf.disable() if useCv else self.ldf.enable()
        self.cvL.enable() if useCv else self.cvL.disable()

    # noinspection PyUnusedLocal,
    def guideCallback(self, *args):
        self.updateGuideGraph()

    def useTheme(self):
        style = ttk.Style(self)
        style.theme_use(self.themeNameVar.get())

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
        for w in self.forceUpdOnThemeWidget:
            w.config(background=fbgc, foreground=fgc, insertbackground=fgc)

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

        grays = [f"gray{i}" for i in [90, 80, 70]] if THEMES[self.themeNameVar.get()] else [f"gray{i}" for i in [16, 23, 30]]
        self.errorText.tag_configure("gun", background=grays[0])
        self.errorText.tag_configure("recoilless", background=grays[0])
        self.errorText.tag_configure("optimize_gun", background=grays[1])
        self.errorText.tag_configure("optimize_recoilless", background=grays[1])
        # self.errorText.tag_configure("minimize", background=grays[2])

        try:
            for fig in (self.fig, self.geomFig, self.auxFig, self.guideFig):
                fig.set_facecolor(bgc)

            for ax in (self.ax, self.axv, self.axP, self.geomAx, self.auxAx, self.auxAxH, self.guideAx):
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
            title=self.getLocStr("exportGraphLabel"),
            filetypes=(("Portable Network Graphics", "*.png"),),
            defaultextension=".png",
            initialfile=filenameize(self.getDescriptive() + " " + save),
        )

        if fileName == "":
            messagebox.showinfo(self.getLocStr("excTitle"), self.getLocStr("cancelMsg"))
            return
        try:
            if save == "main":
                fig = self.fig
            elif save == "aux":
                fig = self.auxFig
            elif save == "guide":
                fig = self.guideFig
            elif save == "geom":
                fig = self.geomFig
            else:
                raise ValueError("unknown save destination")

            with mpl.rc_context(CONTEXT):
                fig.savefig(fileName, transparent=True, dpi=600)

            messagebox.showinfo(self.getLocStr("sucTitle"), self.getLocStr("savedLocMsg") + f" {fileName:}")

        except Exception as e:
            messagebox.showinfo(self.getLocStr("excTitle"), str(e))


# noinspection PyPep8Naming
def calculate(jobQueue, progressQueue, logQueue, kwargs):
    # rootLogger = logging.getLogger()
    rootLogger.addHandler(QueueHandler(logQueue))
    rootLogger.info("calculation started.")

    gunType = kwargs["typ"]
    constrain = kwargs["con"]
    optimize = kwargs["opt"]
    lock = kwargs["lock"]
    # debug = kwargs["deb"]

    try:
        gun, gunResult = None, None
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
        rootLogger.info("calculation concluded successfully.")

    except Exception:
        gun, gunResult = None, None
        exc_type, exc_value, exc_traceback = sys.exc_info()
        rootLogger.error("exception while calculating:")
        rootLogger.error("".join(traceback.format_exception(exc_type, exc_value, exc_traceback)))
    finally:
        # noinspection PyUnboundLocalVariable
        jobQueue.put((kwargs, gun, gunResult))


def guide(guideJobQueue, progressQueue, logQueue, kwargs):
    # rootLogger = logging.getLogger()
    rootLogger.addHandler(QueueHandler(logQueue))
    rootLogger.info("guidance diagram calculation started")

    guideResults = None
    try:
        guideResults = guideGraph(**kwargs)
        rootLogger.info("guidance diagram calculation concluded successfully.")

    except Exception:
        guideResults = None
        exc_type, exc_value, exc_traceback = sys.exc_info()
        rootLogger.error("exception while calculating guidance diagram:")
        rootLogger.error("".join(traceback.format_exception(exc_type, exc_value, exc_traceback)))
    finally:
        # noinspection PyUnboundLocalVariable
        guideJobQueue.put(guideResults)


def main(loc: str = None):
    multiprocessing.freeze_support()
    rootLogger.info("Initializing")

    # this tells windows that our program will handle scaling ourselves
    winRelease = platform.release()
    if winRelease in ("8", "10"):
        windll.shcore.SetProcessDpiAwareness(1)
    elif winRelease in ("7", "Vista"):
        windll.user32.SetProcessDPIAware()
    else:
        print("Unknown release: ", winRelease, ", skipping DPI handling")

    if not loc:
        loc = locale.windows_locale[windll.kernel32.GetUserDefaultUILanguage()]

    # this allows us to set our own taskbar icon
    windll.shell32.SetCurrentProcessExplicitAppUserModelID("Phoenix.Interior.Ballistics.Solver")

    loadfont(resolvepath("ui/sarasa-fixed-sc-regular.ttf"), False, True)
    font_manager.fontManager.addfont(resolvepath("ui/sarasa-fixed-sc-regular.ttf"))

    root = Tk()
    root.option_add("*tearOff", False)
    root.iconbitmap(resolvepath("ui/logo.ico"))

    dpi = root.winfo_fpixels("1i")
    # Tk was originally developed for a dpi of 72
    root.tk.call("tk", "scaling", 1.0 * dpi / 72.0)
    root.tk.call("lappend", "auto_path", resolvepath("ui/tksvg0.14"))
    root.tk.call("lappend", "auto_path", resolvepath("ui/awthemes-10.4.0"))

    root.title("PIBS v" + __version__)
    menubar = Menu(root)
    root.config(menu=menubar)

    root.bind("<Escape>", lambda event: root.state("normal"))
    root.bind("<F11>", lambda event: root.state("zoomed"))
    InteriorBallisticsFrame(root, root, menubar, dpi, defaultLang="English" if loc != "zh_CN" else "中文")
    root.minsize(root.winfo_width(), root.winfo_height())  # set minimum size
    root.mainloop()


if __name__ == "__main__":
    main()
