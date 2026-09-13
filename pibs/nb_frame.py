from __future__ import annotations

import logging
from tkinter import Text, ttk
from tkinter.ttk import Frame, Notebook
from typing import TYPE_CHECKING

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

from . import FONTNAME, FONTSIZE, THEMES
from .ballistics.gun import Gun, GunResult
from .ballistics.prop import Propellant
from .config import SimulationConfig
from .guidegraph import GuideResults
from .localized import LocalizedFrame
from .plot_manager import PlotManager
from .tbl_frame import TableFrame
from .theme import ThemedMixin

if TYPE_CHECKING:
    from .ib_frame import InteriorBallisticsFrame


class NotebookFrame(ThemedMixin, LocalizedFrame):
    def __init__(
        self,
        master: InteriorBallisticsFrame,
        *args,
        font,
        default_lang,
        localization_dict,
        **kwargs,
    ):
        super().__init__(
            master, *args, font=font, default_lang=default_lang, localization_dict=localization_dict, **kwargs
        )

        self.master: InteriorBallisticsFrame = master

        self.font = font

        self.columnconfigure(0, weight=1)
        self.rowconfigure(0, weight=1)
        ## center frame
        self.tab_parent = Notebook(self)

        self.tab_parent.grid(row=0, column=0, sticky="nsew")
        self.tab_parent.columnconfigure(0, weight=1)
        self.tab_parent.rowconfigure(0, weight=1)

        def setup_frame():
            frame = Frame(self.tab_parent)
            frame.grid(row=0, column=0, sticky="nsew")
            frame.rowconfigure(0, weight=1)
            frame.columnconfigure(0, weight=1)
            return frame

        self.desc_tab, self.plot_tab, self.table_tab, self.guide_tab, self.error_tab = (setup_frame() for _ in range(5))

        self.plot_tab.rowconfigure(0, weight=3)
        self.plot_tab.rowconfigure(1, weight=1)

        self.plot_tab.columnconfigure(0, weight=1)
        self.plot_tab.columnconfigure(1, weight=1)

        self.tab_parent.add(self.desc_tab, text=self.get_loc_str("descTab"))
        self.tab_parent.add(self.plot_tab, text=self.get_loc_str("plotTab"))
        self.tab_parent.add(self.table_tab, text=self.get_loc_str("tableTab"))
        self.tab_parent.add(self.guide_tab, text=self.get_loc_str("guideTab"))
        self.tab_parent.add(self.error_tab, text=self.get_loc_str("errorTab"))

        self.tab_parent.enable_traversal()

        ### desc frame
        desc_frm = Frame(self.desc_tab)
        desc_frm.grid(row=0, column=0, sticky="nsew", padx=2, pady=2)
        desc_frm.columnconfigure(0, weight=1)
        desc_frm.rowconfigure(0, weight=1)

        desc_scroll = ttk.Scrollbar(desc_frm, orient="vertical")
        desc_scroll.grid(row=0, column=1, sticky="nsew")
        self.description = Text(
            desc_frm,
            wrap="word",
            height=0,
            width=0,
            yscrollcommand=desc_scroll.set,
            font=(FONTNAME, FONTSIZE),
            undo=True,
            maxundo=-1,
        )
        self.description.grid(row=0, column=0, sticky="nsew")
        desc_scroll.config(command=self.description.yview)
        self.force_update_on_theme_widget.append(self.description)

        ### table frame
        self.table_frame = TableFrame(
            self.table_tab,
            font=self.font,
            default_lang=default_lang,
            localization_dict=localization_dict,
            lang_var=self.lang_var,
        )

        self.table_frame.grid(row=0, column=0, sticky="nsew", padx=2, pady=2)

        ## error frame
        error_frm = self.add_localized_label_frame(self.error_tab, label_loc_key="errFrmLabel")
        error_frm.grid(row=0, column=0, columnspan=3, sticky="nsew", padx=2, pady=2)
        error_frm.columnconfigure(0, weight=1)
        error_frm.rowconfigure(0, weight=1)

        err_scroll = ttk.Scrollbar(error_frm, orient="vertical")
        err_scroll.grid(row=0, column=1, sticky="nsew")
        self.error_text = Text(
            error_frm, yscrollcommand=err_scroll.set, wrap="word", height=0, width=0, font=(FONTNAME, FONTSIZE)
        )
        self.error_text.grid(row=0, column=0, sticky="nsew")

        self.error_text.tag_configure(str(logging.DEBUG), foreground="tan")
        self.error_text.tag_configure(str(logging.WARNING), foreground="orange")
        self.error_text.tag_configure(str(logging.ERROR), foreground="orangered")
        self.error_text.tag_configure(str(logging.CRITICAL), foreground="red")

        err_scroll.config(command=self.error_text.yview)
        self.force_update_on_theme_widget.append(self.error_text)

        ### plot frame
        plot_label_frame = self.add_localized_label_frame(
            self.plot_tab, label_loc_key="plotFrmLabel", tooltip_loc_key="plotText"
        )
        plot_label_frame.grid(row=0, column=0, columnspan=2, sticky="nsew", padx=2, pady=2)

        for i in range(3):
            plot_label_frame.columnconfigure(i, weight=1)

        j = 1
        self.plot_avg_p, self.plot_base_p, self.plot_breech_p = (
            self.add_localized_label_check(
                parent=plot_label_frame, label_loc_key=label, desc_label_key=None, row=j, col=k
            )
            for k, label in enumerate(("plotAvgP", "plotBaseP", "plotBreechP"))
        )

        j += 1
        self.plot_vel, self.plot_nozzle_v, self.plot_burnup = (
            self.add_localized_label_check(
                parent=plot_label_frame, label_loc_key=label, desc_label_key=None, row=j, col=k
            )
            for k, label in enumerate(("plotVel", "plotNozzleV", "plotBurnup"))
        )

        j += 1
        self.plot_stag_p, self.plot_stag_l, self.plot_eta = (
            self.add_localized_label_check(
                parent=plot_label_frame, label_loc_key=label, desc_label_key=None, row=j, col=k
            )
            for k, label in enumerate(("plotStagP", "plotStagL", "plotEta"))
        )

        plot_label_frame.columnconfigure(0, weight=1)
        plot_label_frame.rowconfigure(0, weight=1)

        plot_place_frm = Frame(plot_label_frame)
        plot_place_frm.grid(row=0, column=0, sticky="nsew", columnspan=3)

        fig = Figure(dpi=None, figsize=None, layout="constrained")
        self.fig_canvas = FigureCanvasTkAgg(fig, master=plot_place_frm)
        self.fig_canvas.get_tk_widget().place(relheight=1, relwidth=1)

        aux_plot_label_frame = self.add_localized_label_frame(
            self.plot_tab, label_loc_key="auxFrmLabel", tooltip_loc_key="auxText"
        )
        aux_plot_label_frame.grid(row=1, column=0, sticky="nsew", padx=2, pady=2)

        for i in range(2):
            aux_plot_label_frame.columnconfigure(i, weight=1)

        j = 1
        k = 0
        self.trace_hull, k = (
            self.add_localized_label_check(
                parent=aux_plot_label_frame, row=j, col=k, label_loc_key="traceHull", default=True, desc_label_key=None
            ),
            k + 1,
        )

        self.trace_press, k = (
            self.add_localized_label_check(
                parent=aux_plot_label_frame, row=j, col=k, label_loc_key="tracePress", desc_label_key=None
            ),
            k + 1,
        )

        aux_plot_label_frame.columnconfigure(0, weight=1)
        aux_plot_label_frame.rowconfigure(0, weight=1)

        aux_plot_frame = Frame(aux_plot_label_frame)
        aux_plot_frame.grid(row=0, column=0, sticky="nsew", columnspan=2)

        aux_fig = Figure(dpi=None, layout="constrained")
        self.aux_canvas = FigureCanvasTkAgg(aux_fig, master=aux_plot_frame)
        self.aux_canvas.get_tk_widget().place(relwidth=1, relheight=1)

        ## geom plot
        geom_plot_frm = self.add_localized_label_frame(
            self.plot_tab, label_loc_key="σ(Z)", style="SubLabelFrame.TLabelframe", tooltip_loc_key="geomPlotText"
        )
        geom_plot_frm.grid(row=1, column=1, sticky="nsew", padx=2, pady=2)
        geom_fig = Figure(dpi=None, layout="constrained")
        self.geom_canvas = FigureCanvasTkAgg(geom_fig, master=geom_plot_frm)
        self.geom_canvas.get_tk_widget().place(relheight=1, relwidth=1)

        ## guide plot
        self.guide_tab.rowconfigure(0, weight=1)
        self.guide_tab.columnconfigure(0, weight=1)

        plot_label_frame = self.add_localized_label_frame(self.guide_tab, label_loc_key="guideFrmLabel")
        plot_label_frame.grid(row=0, column=0, columnspan=2, sticky="nsew", padx=2, pady=2)
        guide_fig = Figure(dpi=None, layout="constrained")
        self.guide_canvas = FigureCanvasTkAgg(guide_fig, master=plot_label_frame)
        self.guide_canvas.get_tk_widget().place(relheight=1, relwidth=1)

        plot_label_frame.columnconfigure(0, weight=1)

        guide_plot_option_frame = Frame(self.guide_tab)
        guide_plot_option_frame.grid(row=1, column=0, sticky="nsew")

        self.guide_plot_travel, self.guide_plot_volume, self.guide_plot_burnout, self.guide_chamber_ruler = (
            self.add_localized_label_check(
                guide_plot_option_frame, label_loc_key=label_loc_key, desc_label_key=None, row=i, col=0
            )
            for i, label_loc_key in enumerate(
                ("guidePlotTravel", "guidePlotVolume", "guidePlotBurnout", "guideChamberRuler")
            )
        )
        guide_plot_option_frame.columnconfigure(0, weight=1)

        guide_input_frame = self.add_localized_label_frame(self.guide_tab, label_loc_key="guideCtrlFrmLabel")
        guide_input_frame.grid(row=1, column=1, sticky="nsew", padx=2, pady=2)
        guide_input_frame.columnconfigure(0, weight=1)

        (
            self.guide_step_lf,
            self.guide_min_cmr,
            self.guide_max_cmr,
            self.guide_step_cmr,
        ) = (
            self.add_localized_3_input(
                guide_input_frame,
                label_loc_key=locKey,
                desc_label_key=None,
                default=default,
                unit_text=unit,
                row=j,
                dtype=float,
                label_width=25,
            )
            for j, (locKey, default, unit) in enumerate(
                (
                    ("stepLFLabel", "5.0", "%"),
                    ("minCMRLabel", "0.05", ""),
                    ("maxCMRLabel", "1.00", ""),
                    ("stepCMRLabel", "0.05", ""),
                )
            )
        )

        self.plot_manager = PlotManager(self)

    def use_theme(self):
        super().use_theme()
        grays = (
            [f"gray{i}" for i in [90, 80, 70]]
            if THEMES[self.theme_name_var.get()]["is_light"]
            else [f"gray{i}" for i in [15, 25, 35]]
        )

        self.error_text.tag_configure("base_gun", background=grays[0])
        self.error_text.tag_configure("gun", background=grays[0])
        self.error_text.tag_configure("recoilless", background=grays[0])
        self.error_text.tag_configure("constrained", background=grays[1])
        self.error_text.tag_configure("constrained_gun", background=grays[1])
        self.error_text.tag_configure("constrained_recoilless", background=grays[1])
        self.error_text.tag_configure("guidegraph", background=grays[2])

        self.plot_manager.update_main_plot()
        self.plot_manager.update_aux_plot()
        self.plot_manager.update_guide_graph()
        self.plot_manager.update_geom_plot()

    @property
    def config(self) -> SimulationConfig | None:
        """Simulation configuration from master."""
        return self.master.config

    @property
    def gun(self) -> Gun | None:
        """The gun object from master."""
        return self.master.gun

    @property
    def gun_result(self) -> GunResult | None:
        """Gun simulation result from master."""
        return self.master.gun_result

    @property
    def guide_results(self) -> GuideResults | None:
        """Guide graph result from master."""
        return self.master.guide_results

    @property
    def prop(self) -> Propellant | None:
        return self.master.prop

    def change_lang(self):
        super().change_lang()
        self.tab_parent.tab(self.desc_tab, text=self.get_loc_str("descTab"))
        self.tab_parent.tab(self.plot_tab, text=self.get_loc_str("plotTab"))
        self.tab_parent.tab(self.table_tab, text=self.get_loc_str("tableTab"))
        self.tab_parent.tab(self.error_tab, text=self.get_loc_str("errorTab"))
        self.tab_parent.tab(self.guide_tab, text=self.get_loc_str("guideTab"))

        self.plot_manager.update_main_plot()
        self.plot_manager.update_aux_plot()
        self.plot_manager.update_guide_graph()

        self.table_frame.change_lang()

    def set_description(self, description: str):
        self.description.delete(1.0, "end")
        self.description.insert("end", description)
        self.description.edit_reset()

    def on_state_change(self, type_option: str) -> None:
        self.plot_manager.on_state_change(type_option)
