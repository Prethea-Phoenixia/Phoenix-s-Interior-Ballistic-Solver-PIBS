from __future__ import annotations

import json
import logging
import re
import sys
import traceback
from enum import Enum
from logging.handlers import QueueListener
from math import ceil, log10
from multiprocessing import Process, Queue
from pathlib import Path
from tkinter import IntVar, Menu, StringVar, Text, filedialog, messagebox, ttk
from tkinter.font import Font
from tkinter.ttk import Frame
from typing import Literal

import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

from .notebook_frame import NotebookFrame
from . import THEMES, FONTNAME, BOLDSIZE, FONTSIZE, DESCRIPTION
from . import log_formatter
from . import root_logger
from .ballistics import JSONable
from .ballistics import MIN_BARR_VOLUME, MIN_PROJ_TRAVEL, DOMAIN_LEN, DOMAIN_TIME, CONVENTIONAL, RECOILLESS
from .ballistics import POINT_PEAK_AVG, POINT_PEAK_BREECH, POINT_PEAK_SHOT, POINT_PEAK_STAG
from .ballistics import (
    SOL_LAGRANGE,
    SOL_MAMONTOV,
    SOL_PIDDUCK,
    Composition,
    Geometry,
    Material,
    Propellant,
    SimpleGeometry,
)
from .dispatch import calculate, guide
from .info_frame import InfoFrame
from .localized_widget import LocalizedFrame
from .misc import (
    filenameize,
    format_int_input,
    resolve_path,
    round_sig,
    to_si,
    validate_ce,
    validate_nn,
    validate_pi,
)
from .theme import ThemedMixin
from .tip import create_tool_tip

logger = logging.getLogger(__name__)


class Log(Enum):
    ERROR = logging.ERROR
    WARNING = logging.WARNING


class TextHandler(logging.Handler):
    # This class allows you to log to a Tkinter Text or ScrolledText widget
    # Adapted from Moshe Kaplan: https://gist.github.com/moshekaplan/c425f861de7bbf28ef06

    def __init__(self, text):
        # run the regular Handler __init__
        super().__init__()
        self.setFormatter(log_formatter)
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
        self.text.after_idle(append)


class InteriorBallisticsFrame(ThemedMixin, LocalizedFrame):

    def __init__(self, root, menubar, default_lang, localization_dict, font: Font, os_dark: bool, debug: bool):
        super().__init__(
            root,
            font=font,
            dpi=root.dpi,
            menubar=menubar,
            default_lang=default_lang,
            localization_dict=localization_dict,
            os_dark=os_dark,
        )

        validation_nn = self.register(validate_nn)
        validation_pi = self.register(validate_pi)
        validation_ce = self.register(validate_ce)

        self.font = font
        self.job_queue, self.guide_job_queue, self.log_queue = Queue(), Queue(), Queue()
        self.process, self.guide_process, self.gun_result = None, None, None
        self.root = root

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
        self.data_menu = data_menu
        self.theme_menu = theme_menu
        self.debug_menu = debug_menu

        self.debug = IntVar(value=int(debug))

        design_menu.add_command(label=self.get_loc_str("saveLabel"), command=self.save, accelerator="Ctrl+S")
        self.root.bind("<Control-s>", lambda *_: self.save())
        self.root.bind("<Control-S>", lambda *_: self.save())

        design_menu.add_command(label=self.get_loc_str("loadLabel"), command=self.load_gun, accelerator="Ctrl+L")
        self.root.bind("<Control-l>", lambda *_: self.load_gun())
        self.root.bind("<Control-L>", lambda *_: self.load_gun())

        design_menu.add_command(
            label=self.get_loc_str("loadPresetLabel"),
            command=lambda *_: self.load_gun(initial_dir=resolve_path("examples")),
        )

        design_menu.add_command(label=self.get_loc_str("resetLabel"), command=self.reset, accelerator="Ctrl+N")
        self.root.bind("<Control-N>", lambda *_: self.reset())
        self.root.bind("<Control-n>", lambda *_: self.reset())

        design_menu.add_command(label=self.get_loc_str("calcLabel"), command=self.on_calculate, accelerator="Ctrl+R")
        self.root.bind("<Control-R>", lambda *_: self.on_calculate())
        self.root.bind("<Control-r>", lambda *_: self.on_calculate())

        design_menu.add_command(label=self.get_loc_str("guideLabel"), command=self.on_guide, accelerator="Ctrl+G")
        self.root.bind("<Control-G>", lambda *_: self.on_guide())
        self.root.bind("<Control-g>", lambda *_: self.on_guide())

        data_menu.add_command(label=self.get_loc_str("exportMain"), command=lambda *_: self.export_graph(save="main"))
        data_menu.add_command(label=self.get_loc_str("exportAux"), command=lambda *_: self.export_graph(save="aux"))
        data_menu.add_command(label=self.get_loc_str("exportGeom"), command=lambda *_: self.export_graph(save="geom"))
        data_menu.add_command(label=self.get_loc_str("exportGuide"), command=lambda *_: self.export_graph(save="guide"))
        data_menu.add_command(label=self.get_loc_str("exportLabel"), command=self.export_table)

        data_menu.add_command(label=self.get_loc_str("reloadPropellant"), command=self.load_propellant)

        for theme_name in THEMES.keys():
            theme_menu.add_radiobutton(
                label=theme_name, variable=self.theme_name_var, value=theme_name, command=self.use_theme
            )

        debug_menu.add_checkbutton(label=self.get_loc_str("enableLabel"), variable=self.debug, onvalue=1, offvalue=0)

        self.prop, self.gun, self.guide_result = None, None, None

        self.columnconfigure(1, weight=1)
        self.rowconfigure(1, weight=1)

        ## nameplate
        dummy_frame = Frame(self)
        dummy_frame.grid(row=0, column=0, sticky="nsew", columnspan=2)

        dummy_frame.columnconfigure(0, weight=1)
        dummy_frame.rowconfigure(0, weight=1)

        name_frm = self.add_localized_label_frame(dummy_frame, label_loc_key="nameFrm")
        name_frm.grid(row=0, column=0, sticky="nsew", padx=2, pady=2)

        name_frm.columnconfigure(0, weight=1)
        name_frm.rowconfigure(0, weight=1)
        self.name_var = StringVar(self, value=self.get_loc_str("newDesign"))
        name_plate = ttk.Entry(name_frm, textvariable=self.name_var, justify="left", font=(FONTNAME, BOLDSIZE))
        name_plate.grid(row=0, column=0, sticky="nsew", padx=2, pady=2, columnspan=2)

        self.info_frame = InfoFrame(
            self, font=self.font, default_lang=default_lang, localization_dict=localization_dict
        )
        self.info_frame.grid(row=1, column=0, sticky="nsew")
        self.info_frame.columnconfigure(0, weight=1)
        self.info_frame.rowconfigure(0, weight=1)

        left_frame = Frame(self)
        left_frame.grid(row=0, column=2, rowspan=2, sticky="nsew")
        left_frame.columnconfigure(0, weight=1)
        left_frame.rowconfigure(0, weight=1)

        specs_frame = self.add_localized_label_frame(left_frame, label_loc_key="specFrmLabel")
        specs_frame.grid(row=0, column=0, sticky="nsew", padx=2, pady=2)

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

        specs_frame.rowconfigure(i, weight=1)

        self.grain_frame = self.add_localized_label_frame(specs_frame, label_loc_key="grainFrmLabel")
        self.grain_frame.grid(row=i, column=0, columnspan=3, sticky="nsew", padx=2, pady=2)
        i += 1

        self.grain_frame.columnconfigure(0, weight=1)
        self.grain_frame.rowconfigure(0, weight=1)

        j = 0
        geom_plot_frm = self.add_localized_label_frame(
            self.grain_frame, label_loc_key="σ(Z)", style="SubLabelFrame.TLabelframe", tooltip_loc_key="geomPlotText"
        )
        geom_plot_frm.grid(row=j, column=0, columnspan=3, sticky="nsew", padx=2, pady=2)
        geom_fig = Figure(dpi=None, layout="constrained")
        self.geom_canvas = FigureCanvasTkAgg(geom_fig, master=geom_plot_frm)
        self.geom_canvas.get_tk_widget().place(relheight=1, relwidth=1)

        j += 1

        self.main_geom = self.add_localized_dropdown(
            parent=self.grain_frame, str_obj_dict=Geometry.get_desc_geometry_dict(), desc_label_key="Grain Geometry"
        )
        self.main_geom.grid(row=j, column=0, columnspan=3, sticky="nsew", padx=2, pady=2)

        j += 1

        self.web_mm, j = (
            self.add_localized_3_input(
                parent=self.grain_frame,
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
                parent=self.grain_frame,
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
                parent=self.grain_frame,
                row=j,
                desc_label_key="1/β",
                unit_text="x",
                default="10.0",
                validation=validation_nn,
                dtype=float,
            ),
            j + 1,
        )

        self.use_aux_grain = self.add_localized_label_check(
            parent=self.grain_frame,
            label_loc_key="useAuxGrainLabel",
            desc_label_key="useAuxGrainLabel",
            default=False,
            skip_grid=True,
        )

        self.aux_grain_frm = self.add_localized_label_frame(
            self.grain_frame, labelwidget=self.use_aux_grain.check_widget
        )
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

        self.swap_button = ttk.Button(self.grain_frame, text=self.get_loc_str("swapLabel"), command=self.swap)
        self.swap_button.grid(row=j, column=0, columnspan=3, sticky="nsew", padx=2, pady=2)

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

        self.use_material = self.add_localized_label_check(
            specs_frame,
            label_loc_key="useMaterialLabel",
            desc_label_key="useMaterialLabel",
            default=False,
            skip_grid=True,
        )

        material_frame = self.add_localized_label_frame(specs_frame, labelwidget=self.use_material.check_widget)
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

        self.in_atmos = self.add_localized_label_check(
            parent=left_frame, label_loc_key="atmosLabel", desc_label_key="atmosLabel", skip_grid=True
        )
        environment_frame = self.add_localized_label_frame(left_frame, labelwidget=self.in_atmos.check_widget)
        environment_frame.grid(row=i, column=0, sticky="nsew", padx=2, pady=2)
        i += 1

        environment_frame.columnconfigure(0, weight=1)
        j = 0
        self.amb_p, j = (
            self.add_localized_3_input(
                parent=environment_frame,
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
                parent=environment_frame,
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
                parent=environment_frame,
                row=j,
                label_loc_key="ambGamLabel",
                default="1.400",
                validation=validation_nn,
                dtype=float,
            ),
            j + 1,
        )

        mid_frame = Frame(self)
        mid_frame.grid(row=0, column=3, rowspan=2, sticky="nsew")
        mid_frame.columnconfigure(0, weight=1)

        i = 0
        mid_frame.rowconfigure(i, weight=1)

        ### propellant frame
        prop_frame = self.add_localized_label_frame(mid_frame, label_loc_key="propFrmLabel")
        prop_frame.grid(row=i, column=0, columnspan=3, sticky="nsew", padx=2, pady=2)
        i += 1
        prop_frame.rowconfigure(1, weight=1)
        prop_frame.columnconfigure(0, weight=1)
        j = 0
        self.drop_prop = self.add_localized_dropdown(
            parent=prop_frame,
            str_obj_dict=Composition.read_file(resolve_path("ballistics/resource/propellants.csv")),
            desc_label_key="propFrmLabel",
            tooltip_loc_key="specsText",
        )
        self.drop_prop.grid(row=j, column=0, columnspan=2, sticky="nsew", padx=2, pady=2)
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

        self.use_combustible = self.add_localized_label_check(
            prop_frame,
            label_loc_key="combustibleLabel",
            tooltip_loc_key="combustibleText",
            skip_grid=True,
            default=False,
        )

        combustible_frame = self.add_localized_label_frame(prop_frame, labelwidget=self.use_combustible.check_widget)
        combustible_frame.grid(row=j, column=0, columnspan=2, sticky="nsew", padx=2, pady=2)
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

        sol_frm = self.add_localized_label_frame(mid_frame, label_loc_key="solFrmLabel")
        sol_frm.grid(row=i, column=0, sticky="nsew", padx=2, pady=2)
        i += 1
        sol_frm.columnconfigure(0, weight=1)

        self.drop_gradient = self.add_localized_dropdown(
            parent=sol_frm,
            str_obj_dict={solution: solution for solution in (SOL_LAGRANGE, SOL_PIDDUCK, SOL_MAMONTOV)},
            desc_label_key="solFrmLabel",
        )
        self.drop_gradient.grid(row=0, column=0, columnspan=2, sticky="nsew", padx=2, pady=2)

        op_frm = self.add_localized_label_frame(mid_frame, label_loc_key="opFrmLabel")
        op_frm.grid(row=4, column=0, sticky="nsew", padx=2, pady=2)
        op_frm.columnconfigure(0, weight=1)

        self.use_cons = self.add_localized_label_check(
            parent=op_frm,
            default=False,
            label_loc_key="consButton",
            desc_label_key="consButton",
            tooltip_loc_key="useConsText",
            skip_grid=True,
        )

        cons_frm = self.add_localized_label_frame(op_frm, labelwidget=self.use_cons.check_widget)
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

        self.opt, j = (
            self.add_localized_label_check(
                parent=cons_frm,
                row=j,
                columnspan=3,
                default=False,
                desc_label_key="optButton",
                label_loc_key="optButton",
                tooltip_loc_key="optText",
            ),
            j + 1,
        )

        self.drop_opt_tgt = self.add_localized_dropdown(
            parent=cons_frm,
            str_obj_dict={MIN_BARR_VOLUME: MIN_BARR_VOLUME, MIN_PROJ_TRAVEL: MIN_PROJ_TRAVEL},
            desc_label_key="optTgtLabel",
        )
        self.drop_opt_tgt.grid(row=j, column=0, columnspan=3, sticky="nsew", padx=2, pady=2)
        j += 1

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
                label_loc_key="iniWebLabel",
                unit_text="μm",
                default="100.0",
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

        self.calc_button = ttk.Button(op_frm, text=self.get_loc_str("calcLabel"), command=self.on_calculate)
        self.calc_button.grid(row=i, column=0, columnspan=3, sticky="nsew", padx=2, pady=2)

        op_frm.rowconfigure(i, weight=1)
        i += 1

        self.calc_button_tip = StringVar(value=self.get_loc_str("calcButtonText"))
        create_tool_tip(self.calc_button, self.calc_button_tip, font=self.font)

        self.use_cons.trace_add("write", self.on_state_change)
        self.lock_Lg.trace_add("write", self.on_state_change)
        self.opt.trace_add("write", self.on_state_change)
        self.use_aux_grain.trace_add("write", self.on_state_change)
        self.in_atmos.trace_add("write", self.on_state_change)
        self.type_optn.trace_add("write", self.on_state_change)
        self.use_material.trace_add("write", self.on_state_change)
        self.use_combustible.trace_add("write", self.on_state_change)

        self.main_geom.trace_add("write", self.update_geom)
        self.aux_geom.trace_add("write", self.update_geom)
        self.drop_prop.trace_add("write", self.update_spec)

        for entry in (
            *(self.main_geom, self.aux_geom, self.use_aux_grain, self.drop_prop, self.grain_r1, self.grain_r2),
            *(self.web_mm, self.aux_mass_ratio, self.aux_web_ratio, self.aux_grain_r1, self.aux_grain_r2),
            *(self.use_combustible, self.combustible_force_kJ__kg, self.combustible_mass_kg, self.chg_kg),
            self.force_fudge,
        ):
            entry.trace_add("write", self.propellant_callback)

        self.notebook_frame = NotebookFrame(
            self,
            font=self.font,
            dpi=root.dpi,
            default_lang=default_lang,
            localization_dict=localization_dict,
            lang_var=self.lang_var,
            on_guide_func=self.on_guide,
        )
        self.notebook_frame.grid(row=1, column=1, sticky="nsew", padx=0, pady=0)

        root.protocol("WM_DELETE_WINDOW", self.quit)
        self.use_theme()
        self.reset_entries()  # this fires off traces to initially set the UI into a consistent state.
        self.t_lid = None

        text_handler = TextHandler(self.notebook_frame.error_text)
        root_logger.addHandler(text_handler)
        logger.info("text handler attached to root logger.")

        self.listener = QueueListener(self.log_queue, text_handler)
        self.listener.start()

        logger.info("text handler attached to subprocess log queue listener.")

        self.timed_loop()

    @staticmethod
    def handle_error_wrapper(level: Literal[Log.ERROR, Log.WARNING]):
        def decorator(func):
            def handled_func(self, *args, **kwargs):
                try:
                    return func(self, *args, **kwargs)
                except Exception as e:
                    self.handle_errors(e, level.value)

            return handled_func

        return decorator

    @staticmethod
    def lock_out(func):
        """Manipulating data while a solution or guide graph is being computed can lead to inconsistent program
        state. This decorator locks out execution of a certain program route during said sensitive phases."""

        def handled_func(self, *args, **kwargs):
            if self.process or self.guide_process:
                return None
            else:
                return func(self, *args, **kwargs)

        return handled_func

    def handle_errors(self, exception: Exception, level: int = logging.WARNING):
        if self.debug.get():
            exc_type, exc_value, exc_traceback = sys.exc_info()
            logger.log(level, "".join(traceback.format_exception(exc_type, exc_value, exc_traceback)))
        else:
            logger.log(level, str(exception))

    def _reset_ui(self):
        """Restore buttons and disinhibit widgets after a computation completes."""
        self.calc_button.config(state="normal")
        self.notebook_frame.guide_button.config(state="normal")
        self.swap_button.config(state="normal")
        for loc in self.localized_widgets:
            try:
                loc.disinhibit()
            except AttributeError:
                pass

    def timed_loop(self):
        if self.process:
            self.get_value()

        if self.guide_process:
            self.get_guide()

        self.t_lid = self.root.after(100, self.timed_loop)

    @lock_out
    def reset(self, *_):
        self.reset_entries()

    def quit(self):
        if self.process:
            self.process.terminate()

        if self.guide_process:
            self.guide_process.terminate()

        self.listener.stop()

        if self.t_lid:
            self.root.after_cancel(self.t_lid)

        super().quit()
        self.root.quit()

    def get_normalized_name(self) -> str:
        """
        Return the given string converted to a string that can be used for a clean
        filename. Remove leading and trailing spaces; convert other spaces to
        underscores; and remove anything that is not an alphanumeric, dash,
        underscore, or dot.

        Adapted from Django:
        https://github.com/django/django/blob/main/django/utils/text.py
        """
        s = str(self.name_var.get()).strip().replace(" ", "_")
        s = re.sub(r"(?u)[^-\w.]", "", s)
        if s in {"", ".", ".."}:
            messagebox.showinfo(self.get_loc_str("excTitle"), self.get_loc_str("nameIssue"))
            raise ValueError(self.get_loc_str("nameIssue"))
        return s

    @lock_out
    @handle_error_wrapper(level=Log.WARNING)
    def save(self):
        gun = self.gun
        if gun is None:
            messagebox.showinfo(self.get_loc_str("excTitle"), self.get_loc_str("noDataMsg"))
            return

        file_name = filedialog.asksaveasfilename(
            title=self.get_loc_str("saveLabel"),
            filetypes=(("JSON file", "*.json"),),
            defaultextension=".json",
            initialfile=filenameize(self.get_normalized_name()),
        )

        if not file_name:
            return

        loc_val_dict = {
            loc.get_descriptive(): loc.get()
            for loc in self.localized_widgets
            if (hasattr(loc, "get_descriptive") and loc.get_descriptive())
        }

        kvs = {**loc_val_dict, "Description": self.description.get(1.0, "end").strip("\n")}
        with open(file_name, "w", encoding="utf-8") as file:
            json.dump(kvs, file, indent="\t", ensure_ascii=False, sort_keys=True)

        messagebox.showinfo(self.get_loc_str("sucTitle"), self.get_loc_str("savedLocMsg") + " {:}".format(file_name))

    @lock_out
    @handle_error_wrapper(level=Log.WARNING)
    def load_gun(self, initial_dir: str | None = None):
        file_name = filedialog.askopenfilename(
            title=self.get_loc_str("loadLabel"),
            filetypes=(("JSON File", "*.json"),),
            defaultextension=".json",
            initialdir=initial_dir,
        )

        if not file_name:
            return

        self.reset_entries()
        self.name_var.set(Path(file_name).stem)

        loc_dict = {
            loc.get_descriptive(): loc
            for loc in self.localized_widgets
            if (hasattr(loc, "get_descriptive") and loc.get_descriptive())
        }
        with open(file_name, "r", encoding="utf-8") as file:
            file_dict = json.load(file)

        for key, value in file_dict.items():
            try:  # load design settings (bool -> checkboxes, str -> dropdown menus) first
                if isinstance(value, bool) or isinstance(value, str):
                    loc_dict[key].set(value)

            except (KeyError, ValueError):
                pass

        for key, value in file_dict.items():
            try:  # then load the numeral values (int, float -> numeric entries)
                if isinstance(value, float) or isinstance(value, int):
                    loc_dict[key].set(value)
            except (KeyError, ValueError):
                pass

        if DESCRIPTION in file_dict.keys():  # update description from file.
            self.notebook_frame.set_description(file_dict[DESCRIPTION])

        self.on_calculate()

    @lock_out
    @handle_error_wrapper(Log.WARNING)
    def load_propellant(self):
        file_name = filedialog.askopenfilename(
            title=self.get_loc_str("loadLabel"), filetypes=(("Comma Separated Values File", "*.csv"),)
        )
        self.drop_prop.reset(str_obj_dict=Composition.read_file(file_name))

    def reset_entries(self):
        for loc in self.localized_widgets:
            loc.reset() if hasattr(loc, "get_descriptive") else None

        self.notebook_frame.reset_entries()
        self.name_var.set(self.get_loc_str("newDesign"))
        self.notebook_frame.set_description("")
        self.gun, self.gun_result, self.guide_result = None, None, None

    @lock_out
    @handle_error_wrapper(Log.WARNING)
    def export_table(self):
        self.notebook_frame.table_frame.export_table(
            gun_result=self.gun_result,
            acc_exp=int(self.acc_exp.get()),
            normalized_filename=filenameize(self.get_normalized_name()),
        )

    @lock_out
    def change_lang(self):
        super().change_lang()
        self.info_frame.change_lang()

        self.notebook_frame.change_lang()

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
        self.calc_button.config(text=self.get_loc_str("calcLabel"))
        self.swap_button.config(text=self.get_loc_str("swapLabel"))
        self.update_geom_plot()
        self.update_table()

    @property
    def kwargs(self) -> dict:
        kwargs = {}
        try:
            constrain = bool(self.use_cons.get())
            lock = bool(self.lock_Lg.get())
            optimize = bool(self.opt.get())
            debug = bool(self.debug.get())
            atmosphere = bool(self.in_atmos.get())
            autofrettage = bool(self.material_is_af.get())
            material = bool(self.use_material.get())
            gun_type = self.type_optn.get_obj()

            if self.prop is None:
                raise ValueError("Invalid propellant.")

            chamber_volume = float(self.cv_L.get()) * 1e-3

            chambrage = float(self.clr.get())
            charge_mass = float(self.chg_kg.get())
            caliber = float(self.cal_mm.get()) * 1e-3

            gun_length = float(self.tbl_mm.get()) * 1e-3
            load_fraction = charge_mass / chamber_volume / self.prop.rho_p

            kwargs = {
                "opt": optimize,
                "con": constrain,
                "deb": debug,
                "lock": lock,
                "typ": gun_type,
                "dom": self.drop_domain.get_obj(),
                "sol": self.drop_gradient.get_obj(),
                "control": self.p_control.get_obj(),
                "structural_material": (
                    Material(
                        density=float(self.material_density.get()),
                        yield_strength=float(self.material_yield.get()) * 1e6,
                    )
                    if material
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
                "min_cmr": float(self.notebook_frame.guide_min_cmr.get()),
                "max_cmr": float(self.notebook_frame.guide_max_cmr.get()),
                "step_cmr": float(self.notebook_frame.guide_step_cmr.get()),
                "step_lf": float(self.notebook_frame.guide_step_lf.get()) * 1e-2,
                "max_iteration": int(self.max_iter.get()),
                "opt_target": self.drop_opt_tgt.get_obj(),
            }

            if atmosphere:
                kwargs.update(
                    {
                        "ambient_pressure": float(self.amb_p.get()) * 1e3,
                        "ambient_density": float(self.amb_rho.get()),
                        "ambient_adb_index": float(self.amb_gamma.get()),
                    }
                )
            else:
                kwargs.update({"ambient_pressure": 0.0, "ambient_density": 0.0, "ambient_adb_index": 1.0})
        except ValueError as e:
            self.handle_errors(e, level=logging.ERROR)

        return kwargs

    def kwargs_to_json(self) -> dict:
        serialized_kwargs = {}
        for key, value in self.kwargs.items():
            if isinstance(value, JSONable):
                value = json.loads(value.to_json())
            serialized_kwargs[key] = value

        return serialized_kwargs

    def _inhibit_widgets(self):
        for loc in self.localized_widgets:
            try:
                loc.inhibit()
            except AttributeError:
                pass

        self.calc_button.config(state="disabled")
        self.notebook_frame.guide_button.config(state="disabled")
        self.swap_button.config(state="disabled")

    def on_guide(self):
        self.focus()  # remove focus to force widget entry validation

        self.guide_process = Process(target=guide, args=(self.guide_job_queue, self.log_queue, self.kwargs))
        self.guide_process.start()

        self._inhibit_widgets()

    def on_calculate(self):
        self.focus()  # remove focus to force widget entry validation

        self.process = Process(target=calculate, args=(self.job_queue, self.log_queue, self.kwargs))
        self.process.start()

        self._inhibit_widgets()

    def get_value(self):
        kwargs: dict[str, int | float] | None = None
        while not self.job_queue.empty():
            kwargs, self.gun, self.gun_result = self.job_queue.get_nowait()

        if kwargs is None:
            return

        try:
            constrain = kwargs["con"]
            lock = kwargs["lock"]
            optimize = kwargs["opt"]
            sigfig = int(-log10(kwargs["tol"])) + 1
            if self.gun:
                if constrain:
                    web_mm = round_sig(kwargs["web"] * 1e3, n=sigfig)
                    self.web_mm.set(web_mm)
                    if not lock:
                        lg_mm = round_sig(kwargs["length_gun"] * 1e3, n=sigfig)
                        self.tbl_mm.set(lg_mm)
                    if optimize:
                        self.cv_L.set(round_sig(kwargs["chamber_volume"] * 1e3, n=sigfig))
            else:
                pass

        except Exception as e:
            self.handle_errors(e, level=logging.WARNING)

        self.update_stats()
        self.update_table()

        self.notebook_frame.update_fig_plot()
        self.notebook_frame.update_aux_plot()
        self.notebook_frame.update_guide_graph()

        self._reset_ui()
        self.process = None

    def get_guide(self):

        initial_result = self.guide_result
        while not self.guide_job_queue.empty():
            self.guide_result = self.guide_job_queue.get_nowait()
        if initial_result == self.guide_result:
            return

        self.notebook_frame.update_guide_graph()
        self._reset_ui()
        self.guide_process = None

    @handle_error_wrapper(Log.WARNING)
    def update_table(self):
        self.notebook_frame.table_frame.update_table(gun_result=self.gun_result, acc_exp=int(self.acc_exp.get()))

    @handle_error_wrapper(Log.WARNING)
    def update_stats(self):
        self.info_frame.update_stats(self.gun, self.gun_result, int(self.acc_exp.get()))

    def update_spec(self, *_):
        self.specs.config(state="normal")
        compo: Composition = self.drop_prop.get_obj()
        self.specs.delete("1.0", "end")

        if compo.temp_v:
            self.specs.insert(
                "end",
                "{:}: {:>4.0f} K {:}\n".format(
                    self.get_loc_str("TvDesc"), compo.temp_v, self.get_loc_str("isochorDesc")
                ),
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
        for p in (1e6, 10e6, 100e6, 1000e6):
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
                r2.localize("ltarcLabel", "ltarcText")

    def update_geom_plot(self):
        with plt.rc_context(self.context):
            self.geom_canvas.figure.clear()
            self.geom_canvas.figure.set_facecolor(self.context["figure.facecolor"])
            geom_ax = self.geom_canvas.figure.add_subplot(111)

            n = 100

            prop = self.prop
            if prop is not None:
                zb = prop.z_b
                xs = [(i / n) * zb for i in range(n + 1)]
                ys = [prop.f_sigma_z(x) for x in xs]

                xs.append(zb)
                ys.append(prop.f_sigma_z(zb))

                xs.append(xs[-1])
                ys.append(0)

                geom_ax.plot(xs, ys)
                geom_ax.grid(which="major", color="grey", linestyle="dotted")
                geom_ax.minorticks_on()
                geom_ax.set_xlim(left=0, right=prop.z_b)
                geom_ax.xaxis.set_ticks([i * 0.5 for i in range(ceil(min(prop.z_b, 2) / 0.5) + 1)])
                geom_ax.set_ylim(bottom=0, top=max(ys))
                geom_ax.yaxis.set_ticks([i * 0.25 for i in range(ceil(max(ys) / 0.25) + 1)])

            self.geom_canvas.draw_idle()

    @handle_error_wrapper(Log.WARNING)
    def propellant_callback(self, *_):
        """
        updates the propellant object on write to the ratio entry fields
        and, on changing the propellant or geometrical specification.
        """

        try:
            self.prop = Propellant(
                composition=self.drop_prop.get_obj(),
                main_geom=self.main_geom.get_obj(),
                main_r1=self.grain_r1.get(),
                main_r2=self.grain_r2.get(),
                aux_geom=self.aux_geom.get_obj(),
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
            self.update_geom_plot()
        except Exception as e:
            self.prop = None
            self.update_geom_plot()

            raise e

    def on_state_change(self, *_):
        self.notebook_frame.on_state_change(type_option=self.type_optn.get_obj())
        if self.type_optn.get_obj() == CONVENTIONAL:
            self.drop_gradient.enable()

            self.nozz_exp.remove()
            self.nozz_eff.remove()

            self.p_control.reset(
                {point: point for point in (POINT_PEAK_AVG, POINT_PEAK_SHOT, POINT_PEAK_BREECH)}, overwrite=False
            )

        elif self.type_optn.get_obj() == RECOILLESS:
            self.drop_gradient.set_by_obj(SOL_LAGRANGE)
            self.drop_gradient.disable()

            self.nozz_exp.restore()
            self.nozz_eff.restore()

            self.p_control.reset(
                {point: point for point in (POINT_PEAK_AVG, POINT_PEAK_SHOT, POINT_PEAK_STAG, POINT_PEAK_BREECH)},
                overwrite=False,
            )

        if not self.use_cons.get():
            self.v_tgt.disable()
            self.p_tgt.disable()

            self.opt.disable()
            self.lock_Lg.disable()
            self.min_web.disable()
            self.lg_max.disable()
            self.p_control.disable()
            self.tbl_mm.enable()

            self.drop_opt_tgt.disable()

        else:
            if self.lock_Lg.get():
                self.v_tgt.disable()
                self.opt.disable()
                self.tbl_mm.enable()
            else:
                self.v_tgt.enable()
                self.opt.enable()
                self.tbl_mm.disable()

            if self.opt.get():
                self.lock_Lg.disable()
                self.drop_opt_tgt.enable()
            else:
                self.lock_Lg.enable()
                self.drop_opt_tgt.disable()

            self.p_tgt.enable()
            self.min_web.enable()
            self.lg_max.enable()
            self.p_control.enable()

        for entry in (self.aux_grain_r1, self.aux_grain_r2, self.aux_web_ratio, self.aux_mass_ratio, self.aux_geom):
            if self.use_aux_grain.get():
                entry.enable()
            else:
                entry.disable()

        for entry in (self.material_yield, self.material_density, self.material_ssf, self.material_is_af):
            if self.use_material.get():
                entry.enable()
            else:
                entry.disable()

        for entry in (self.amb_p, self.amb_rho, self.amb_gamma):
            if self.in_atmos.get():
                entry.enable()
            else:
                entry.disable()

        for entry in (self.combustible_mass_kg, self.combustible_force_kJ__kg):
            if self.use_combustible.get():
                entry.enable()
            else:
                entry.disable()

    @lock_out
    def use_theme(self):
        super().use_theme()
        self.notebook_frame.use_theme()
        self.update_geom_plot()

    @handle_error_wrapper(level=Log.WARNING)
    def export_graph(self, save: Literal["main", "aux", "guide", "geom"]):
        file_name = filedialog.asksaveasfilename(
            title=self.get_loc_str("exportGraphLabel"),
            filetypes=(("Portable Network Graphics", "*.png"),),
            defaultextension=".png",
            initialfile=filenameize(f"{self.get_normalized_name()}_{save}"),
        )

        if save == "main":
            fig = self.notebook_frame.fig_canvas.figure
        elif save == "aux":
            fig = self.notebook_frame.aux_canvas.figure
        elif save == "guide":
            fig = self.notebook_frame.guide_canvas.figure
        elif save == "geom":
            fig = self.geom_canvas.figure
        else:
            raise ValueError("unknown save target.")

        with plt.rc_context(self.context):
            fig.savefig(file_name, transparent=True, dpi=600)

        messagebox.showinfo(self.get_loc_str("sucTitle"), self.get_loc_str("savedLocMsg") + f" {file_name:}")

    @handle_error_wrapper(level=Log.WARNING)
    def swap(self):
        ## this swaps the primary and auxiliary charge.
        sigfig = int(self.acc_exp.get()) + 1
        # cache the values for swapping.
        main_geom, aux_geom = self.main_geom.get_obj(), self.aux_geom.get_obj()
        main_r1, main_r2 = self.grain_r1.get(), self.grain_r2.get()
        aux_r1, aux_r2 = self.aux_grain_r1.get(), self.aux_grain_r2.get()
        web_ratio, mass_ratio = self.aux_web_ratio.get(), self.aux_mass_ratio.get()

        self.grain_r1.set(aux_r1)
        self.grain_r2.set(aux_r2)
        self.aux_grain_r1.set(main_r1)
        self.aux_grain_r2.set(main_r2)

        self.aux_web_ratio.set(round_sig(1.0 / web_ratio, n=sigfig))
        self.aux_mass_ratio.set(round_sig(1.0 / mass_ratio, n=sigfig))

        self.main_geom.set_by_obj(aux_geom)
        self.aux_geom.set_by_obj(main_geom)
