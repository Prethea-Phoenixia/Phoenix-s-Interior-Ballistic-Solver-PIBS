from __future__ import annotations

import logging
import math
from tkinter import ttk, Text
from tkinter.ttk import Frame, Notebook
from typing import Protocol, runtime_checkable

import matplotlib as mpl
from labellines import labelLines
from matplotlib import pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

from ballistics import Propellant
from .ballistics.recoilless import RecoillessTableEntry
from . import FONTSIZE, FONTNAME, THEMES
from .ballistics import DOMAIN_TIME, DOMAIN_LEN, CONVENTIONAL, RECOILLESS
from .ballistics.gun import GunResult, Gun
from .localized_widget import LocalizedFrame
from .misc import validate_nn, validate_ce
from .table_frame import TableFrame
from .theme import ThemedMixin


@runtime_checkable
class MasterInterface(Protocol):
    """Protocol defining the expected interface for NotebookFrame's master."""

    gun: Gun | None
    """The gun object or None if not computed."""

    kwargs: dict[str, object]
    """Dictionary containing simulation parameters including:
        - design_velocity: float
        - design_pressure: float
        - typ: str (CONVENTIONAL or RECOILLESS)
        - dom: str (DOMAIN_TIME or DOMAIN_LEN)
        - control: str
        - charge_mass: float
        - chamber_volume: float
    """

    gun_result: GunResult | None
    """The gun simulation result, None on startup or invalid compute output."""

    guide_result: list[tuple] | None
    """The guide graph result, None on startup or invalid compute output."""

    prop: Propellant | None
    """The propellant, None if invalid propellant entry state"""


class NotebookFrame(ThemedMixin, LocalizedFrame):
    def __init__(
        self,
        master: MasterInterface | None = None,
        *args,
        font,
        default_lang,
        localization_dict,
        on_guide_func,
        **kwargs,
    ):
        super().__init__(master, *args, default_lang=default_lang, localization_dict=localization_dict, **kwargs)

        if not isinstance(master, MasterInterface):
            raise ValueError
        self.master: MasterInterface = master  # type: ignore[assignment]

        validation_nn = self.register(validate_nn)
        validation_ce = self.register(validate_ce)
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
            height=40,
            width=80,
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
            error_frm, yscrollcommand=err_scroll.set, wrap="word", height=40, width=80, font=(FONTNAME, FONTSIZE)
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

        fig = Figure(dpi=None, layout="constrained")
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
                # desc_label_key=None,
                default=default,
                unit_text=unit,
                validation=validation,
                row=j,
                dtype=float,
            )
            for j, (locKey, default, unit, validation) in enumerate(
                (
                    ("stepLFLabel", "1.0", "%", validation_ce),
                    ("minCMRLabel", "0.05", "", validation_nn),
                    ("maxCMRLabel", "1.00", "", validation_nn),
                    ("stepCMRLabel", "0.05", "", validation_nn),
                )
            )
        )
        self.guide_button = ttk.Button(guide_input_frame, text=self.get_loc_str("guideLabel"), command=on_guide_func)
        self.guide_button.grid(row=j + 1, column=0, columnspan=3, sticky="nsew", padx=2, pady=2)

        for check in (
            *(self.plot_avg_p, self.plot_base_p, self.plot_breech_p, self.plot_stag_p, self.plot_stag_l, self.plot_vel),
            *(self.plot_nozzle_v, self.plot_burnup, self.plot_eta),
        ):
            check.trace_add("write", lambda *_: self.update_fig_plot())

        for check in (self.trace_hull, self.trace_press):
            check.trace_add("write", lambda *_: self.update_aux_plot())

        for entry in (self.guide_plot_travel, self.guide_plot_volume, self.guide_plot_burnout):
            entry.trace_add("write", lambda *_: self.update_guide_graph())
        self.guide_chamber_ruler.trace_add("write", lambda *_: self.update_guide_graph())

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

        self.update_fig_plot()
        self.update_aux_plot()
        self.update_guide_graph()
        self.update_geom_plot()

    @property
    def gun(self) -> Gun | None:
        """The gun object from master."""
        return self.master.gun

    @property
    def kwargs(self) -> dict[str, object]:
        """Simulation parameters dictionary from master."""
        return self.master.kwargs

    @property
    def gun_result(self) -> GunResult | None:
        """Gun simulation result from master."""
        return self.master.gun_result

    @property
    def guide_result(self) -> list[tuple] | None:
        """Guide graph result from master."""
        return self.master.guide_result

    @property
    def prop(self) -> Propellant | None:
        return self.master.prop

    def change_lang(self):
        super().change_lang()
        self.guide_button.config(text=self.get_loc_str("guideLabel"))
        self.tab_parent.tab(self.desc_tab, text=self.get_loc_str("descTab"))
        self.tab_parent.tab(self.plot_tab, text=self.get_loc_str("plotTab"))
        self.tab_parent.tab(self.table_tab, text=self.get_loc_str("tableTab"))
        self.tab_parent.tab(self.error_tab, text=self.get_loc_str("errorTab"))
        self.tab_parent.tab(self.guide_tab, text=self.get_loc_str("guideTab"))

        self.update_fig_plot()
        self.update_aux_plot()
        self.update_guide_graph()

        self.table_frame.change_lang()

    def update_fig_plot(self):

        with plt.rc_context(self.context):
            self.fig_canvas.figure.clear()
            self.fig_canvas.figure.set_facecolor(self.context["figure.facecolor"])
            ax = self.fig_canvas.figure.add_subplot(111)
            ax_p, ax_v = ax.twinx(), ax.twinx()

            if self.gun:
                v_tgt = self.kwargs["design_velocity"]
                p_tgt = self.kwargs["design_pressure"]
                gun_type = self.kwargs["typ"]
                dom = self.kwargs["dom"]

                xs, vs, vxs, pas, pss, pbs, p0s, psis, etas, stags = [], [], [], [], [], [], [], [], [], []

                for entry in self.gun_result.table_data:
                    tag = entry.tag
                    time = entry.time
                    travel = entry.travel
                    psi = entry.burnup
                    v = entry.velocity
                    pb = entry.breech_pressure

                    if isinstance(entry, RecoillessTableEntry):
                        vx = entry.outflow_velocity
                        p0 = entry.stag_pressure
                        eta = entry.outflow_fraction
                        stag = entry.rel_stag_point

                    else:
                        vx, p0, eta, stag = 0, 0, 0, 0

                    p = entry.avg_pressure
                    ps = entry.shot_pressure

                    if tag == self.kwargs["control"]:
                        x_peak = (time * 1e3) if dom == DOMAIN_TIME else travel
                        # noinspection PyTypeChecker
                        ax_p.spines.right.set_position(("data", x_peak))

                    if dom == DOMAIN_TIME:
                        xs.append(time * 1000)
                    elif dom == DOMAIN_LEN:
                        xs.append(travel)

                    vs.append(v)
                    vxs.append(vx)
                    pas.append(p * 1e-6)
                    pss.append(ps * 1e-6)
                    pbs.append(pb * 1e-6)
                    p0s.append(p0 * 1e-6)
                    psis.append(psi)
                    etas.append(eta)
                    stags.append(stag)

                if self.plot_breech_p.get():
                    ax_p.plot(
                        xs,
                        pbs,
                        c="xkcd:goldenrod",
                        label=(self.get_loc_str("figBreech" if gun_type == CONVENTIONAL else "figNozzleP")),
                    )

                if gun_type == RECOILLESS:
                    if self.plot_stag_p.get():
                        ax_p.plot(xs, p0s, "seagreen", label=self.get_loc_str("figStagnation"))

                    if self.plot_stag_l.get():
                        ax.plot(xs, stags, "mediumorchid", label=self.get_loc_str("figStagL"))

                    if self.plot_nozzle_v.get():
                        ax_v.plot(xs, vxs, "steelblue", label=self.get_loc_str("figNozzleV"))

                    if self.plot_eta.get():
                        ax.plot(xs, etas, "maroon", label=self.get_loc_str("figOutflow"))

                if self.plot_avg_p.get():
                    ax_p.plot(xs, pas, "tab:green", label=self.get_loc_str("figAvgP"))

                if self.plot_base_p.get():
                    ax_p.plot(xs, pss, "yellowgreen", label=self.get_loc_str("figShotBase"))

                if gun_type == CONVENTIONAL or gun_type == RECOILLESS:
                    ax_p.axhline(float(p_tgt), c="tab:green", linestyle=":", label=self.get_loc_str("figTgtP"))

                if self.plot_vel.get():
                    ax_v.plot(xs, vs, "tab:blue", label=self.get_loc_str("figShotVel"))
                ax_v.axhline(v_tgt, c="tab:blue", linestyle=":", label=self.get_loc_str("figTgtV"))

                if self.plot_burnup.get():
                    ax.plot(xs, psis, c="crimson", label=self.get_loc_str("figPsi"))

                lines_labeled = []
                for lines, xvals in zip(
                    (ax_p.get_lines(), ax.get_lines(), ax_v.get_lines()),
                    (
                        (0.2 * xs[-1] + 0.8 * x_peak, xs[-1]),
                        (0, xs[-1]),
                        (x_peak, 0.2 * xs[-1] + 0.8 * x_peak),
                        (0, x_peak),
                    ),
                ):
                    labelLines(lines, align=True, xvals=xvals)
                    lines_labeled.append(lines)

                ax.set_xlim(left=0, right=xs[-1])
                pmax = max(pas + pbs + pss + p0s)
                ax_p.set(ylim=(0, pmax * 1.1))
                ax_v.set(ylim=(0, max(vs + vxs) * 1.15))
                ax.set_ylim(bottom=0, top=1.05)

                ax_p.yaxis.set_ticks([v for v in ax_p.get_yticks() if v <= pmax][1:])

                ax.yaxis.tick_right()
                ax_p.yaxis.tick_right()
                ax_v.yaxis.tick_left()

                ax.tick_params(axis="y", colors="tab:red")
                ax_v.tick_params(axis="y", colors="tab:blue")
                ax_p.tick_params(axis="y", colors="tab:green")
                ax.tick_params(axis="x")

                ax_p.yaxis.set_label_position("right")
                ax_p.set_ylabel("MPa")
                ax_p.yaxis.label.set_color("tab:green")

                ax_v.yaxis.set_label_position("left")
                ax_v.set_ylabel("m/s")
                ax_v.yaxis.label.set_color("tab:blue")

                if dom == DOMAIN_TIME:
                    ax.set_xlabel(self.get_loc_str("figTimeDomain"))
                elif dom == DOMAIN_LEN:
                    ax.set_xlabel(self.get_loc_str("figLenDomain"))

            else:
                pass
            self.fig_canvas.draw_idle()

    def update_aux_plot(self):

        with plt.rc_context(self.context):
            self.aux_canvas.figure.clear()
            self.aux_canvas.figure.set_facecolor(self.context["figure.facecolor"])
            aux_ax = self.aux_canvas.figure.add_subplot(111)
            aux_ax_h = aux_ax.twinx()
            aux_ax_h.yaxis.tick_right()

            if self.gun_result:
                p_trace = self.gun_result.pressure_trace
                cmap = mpl.colormaps[THEMES[self.theme_name_var.get()]["cmap"]]
                x_max, y_max, t_min, t_max = 0.0, 0.0, math.inf, 0.0
                for trace in p_trace:

                    if not trace.temperature:
                        continue
                    if trace.temperature > t_max:
                        t_max = trace.temperature
                    elif trace.temperature < t_min:
                        t_min = trace.temperature

                for trace in p_trace[::-1]:

                    tag, t, trace = trace.tag, trace.temperature, trace.pressure_trace

                    if t:
                        v = (t - t_min) / (t_max - t_min)
                        color = cmap(v)
                    else:
                        color = cmap(0.5)
                    linestyle = None
                    alpha = None

                    x, y = zip(*[(ppp.x, ppp.p) for ppp in trace])
                    y = [v * 1e-6 for v in y]
                    x_max = max(x_max, float(max(x)))
                    y_max = max(y_max, max(y))

                    if self.trace_press.get():
                        aux_ax.plot(x, y, c=color, alpha=alpha, ls=linestyle)

                aux_ax.set_xlim(left=0, right=x_max)
                aux_ax.set_ylim(bottom=0, top=y_max * 1.15)

                aux_ax.tick_params(axis="y", colors="tab:green")
                aux_ax_h.tick_params(axis="y", colors="tab:blue")
                aux_ax.tick_params(axis="x")

                aux_ax.set_xlabel(self.get_loc_str("figAuxDomain"))

                aux_ax.yaxis.label.set_color("tab:green")
                aux_ax.set_ylabel("MPa")

                aux_ax_h.yaxis.set_ticks_position("right")
                aux_ax_h.yaxis.set_label_position("right")

                aux_ax_h.yaxis.label.set_color("tab:blue")
                aux_ax_h.set_ylabel("mm")

                h_trace = self.gun_result.outline

                if h_trace is not None and self.trace_hull.get():

                    x_hull = list(entry.x for entry in h_trace)
                    r_in = list(entry.r_in * 1e3 for entry in h_trace)
                    r_out = list(entry.r_ex * 1e3 for entry in h_trace)
                    r_pej = list(entry.r_pej * 1e3 for entry in h_trace)

                    aux_ax_h.fill_between(
                        x_hull, r_in, r_pej, alpha=0.5 if self.trace_press.get() else 0.8, color="tab:orange"
                    )
                    aux_ax_h.fill_between(
                        x_hull, r_pej, r_out, alpha=0.5 if self.trace_press.get() else 0.8, color="tab:blue"
                    )

                    aux_ax.set_xlim(left=min(x_hull))

                aux_ax_h.set_ylim(bottom=0)

            self.aux_canvas.draw_idle()

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
                geom_ax.xaxis.set_ticks([i * 0.5 for i in range(math.ceil(min(prop.z_b, 2) / 0.5) + 1)])
                geom_ax.set_ylim(bottom=0, top=max(ys))
                geom_ax.yaxis.set_ticks([i * 0.25 for i in range(math.ceil(max(ys) / 0.25) + 1)])

            self.geom_canvas.draw_idle()

    def update_guide_graph(self):
        style = ttk.Style(self)
        # bgc = str(style.lookup("TFrame", "background"))
        fgc = str(style.lookup("TFrame", "foreground"))
        cmap = mpl.colormaps[THEMES[self.theme_name_var.get()]["cmap"]].reversed()

        def _scale_bounds(values: list[float]):
            min_val = min([value for value in values if value])
            max_val = max([(value if value else min_val) for value in values])
            min_log = int(math.floor(math.log10(min_val)))
            max_log = int(math.ceil(math.log10(max_val)))
            return min_val, max_val, min_log, max_log

        def get_adaptive_scale(values: list[float]) -> list[float]:
            min_val, _max_val, min_log, max_log = _scale_bounds(values)
            base = math.floor(min_val / 10 ** (min_log - 1)) * 10 ** (min_log - 1)
            levels = []

            for j in range(11):
                levels.append(base + j * 10 ** (min_log - 1))

            for i in range(max_log - min_log):
                for j in range(9):
                    v = 10 ** (min_log + i) + j * 10 ** (min_log + i)
                    if v > levels[-1]:
                        levels.append(v)

            levels.append(10**max_log)
            return levels

        def get_decimal_scale(values: list[float]) -> list[float]:
            min_val, _max_val, min_log, max_log = _scale_bounds(values)
            levels = []

            for i in range(max_log - min_log):
                for j in range(9):
                    v = 10 ** (min_log + i) + j * 10 ** (min_log + i)
                    levels.append(v)

            levels.append(10**max_log)
            return levels

        def get_01_scale(_):
            return [0.1 * (i + 1) for i in range(9)]

        with plt.rc_context(self.context):
            self.guide_canvas.figure.clear()
            self.guide_canvas.figure.set_facecolor(self.context["figure.facecolor"])
            guide_ax = self.guide_canvas.figure.add_subplot(111)

            if self.guide_result and any(line for line in self.guide_result):
                load_densities = list(line[0] for line in self.guide_result)
                charge_masses = list(line[1] for line in self.guide_result)

                delta_max = max(load_densities)
                delta_min = min(load_densities)
                w_max = max(charge_masses)
                w_min = min(charge_masses)
                max_cv = w_max / delta_min
                min_cv = w_min / delta_max

                left_diagonal_chamber_volume = w_min / delta_min
                right_diagonal_chamber_volume = w_max / delta_max

                if self.gun:
                    guide_ax.scatter(
                        float(self.kwargs["charge_mass"]) / float(self.kwargs["chamber_volume"]),
                        self.kwargs["charge_mass"],
                        c=fgc,
                        marker="x",
                        s=FONTSIZE * 4,
                    )

                if self.guide_chamber_ruler.get():

                    result_levels = get_decimal_scale([min_cv, max_cv])

                    for cv in result_levels:
                        if cv < min_cv or cv > max_cv:
                            continue
                        if cv < left_diagonal_chamber_volume:
                            lp = (w_min / cv, w_min)
                        else:
                            lp = (delta_min, cv * delta_min)

                        if cv < right_diagonal_chamber_volume:
                            rp = (delta_max, cv * delta_max)

                        else:
                            rp = (w_max / cv, w_max)

                        guide_ax.plot(*zip(lp, rp), color=fgc, label=f"{cv*1e3:.3g} L")
                        labelLines(guide_ax.get_lines(), drop_label=True)

                guide_ax.set_xlabel(self.get_loc_str("guideLDDomain"))
                guide_ax.set_ylabel(self.get_loc_str("guideCMDomain"))

                titles = []

                for index, show, scaling, levels_func, title_loc_str, linestyle, unit in zip(
                    (3, 4, 5),
                    (self.guide_plot_travel.get(), self.guide_plot_volume.get(), self.guide_plot_burnout.get()),
                    (1, 1000, 1),
                    (get_adaptive_scale, get_adaptive_scale, get_01_scale),
                    ("guideTravelTitle", "guideBVTitle", "guideBurnoutTitle"),
                    ("--", "-.", ":"),
                    (" m", " L", ""),
                ):
                    if not show:
                        continue

                    titles.append(self.get_loc_str(title_loc_str))
                    results = list(line[index] * scaling for line in self.guide_result)
                    result_levels = levels_func(results)
                    # noinspection PyTypeChecker
                    cs = guide_ax.tricontour(
                        load_densities,
                        charge_masses,
                        results,
                        levels=result_levels,
                        cmap=cmap,
                        linestyles=linestyle,
                        vmin=max(min(results), min(result_levels)),
                        vmax=min(max(results), max(result_levels)),
                    )
                    guide_ax.clabel(cs, cs.levels, fontsize=FONTSIZE, fmt=lambda v: f"{v:.4g}{unit}")

                guide_ax.set_title("\n".join(titles))
                guide_ax.set_xlim(min(load_densities), max(load_densities))
                guide_ax.set_ylim(min(charge_masses), max(charge_masses))

            self.guide_canvas.draw_idle()

    def set_description(self, description: str):
        self.description.delete(1.0, "end")
        self.description.insert("end", description)
        self.description.edit_reset()

    def reset_entries(self):
        for loc in self.localized_widgets:
            loc.reset() if hasattr(loc, "get_descriptive") else None

    def on_state_change(self, type_option: str) -> None:

        if type_option == CONVENTIONAL:

            self.plot_nozzle_v.remove()
            self.plot_breech_p.localize("plotBreechP")
            self.plot_stag_p.remove()
            self.plot_stag_l.remove()
            self.plot_eta.remove()

        elif type_option == RECOILLESS:

            self.plot_nozzle_v.restore()
            self.plot_breech_p.localize("plotNozzleP")
            self.plot_stag_p.restore()
            self.plot_stag_p.localize("plotStagP")
            self.plot_stag_l.restore()
            self.plot_stag_l.localize("plotStagL")
            self.plot_eta.restore()
            self.plot_eta.localize("plotEtaEsc")

        else:
            raise ValueError
