from __future__ import annotations

import math
from tkinter import ttk
from typing import TYPE_CHECKING

import matplotlib as mpl
from labellines import labelLines
from matplotlib import pyplot as plt

from . import BOLDSIZE, FONTNAME, FONTSIZE, THEMES, FigureType, Theme
from .ballistics import Domain, GunType
from .ballistics.recoilless import RecoillessTableEntry

if TYPE_CHECKING:
    from .ib_frame import InteriorBallisticsFrame
    from .nb_frame import NotebookFrame


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


class PlotManager:
    """Manages all matplotlib plot updates."""

    def __init__(self, frame: NotebookFrame):
        self.frame = frame
        self._setup_traces()

    def _setup_traces(self):
        """Set up trace bindings for plot checkboxes."""
        f = self.frame

        # Main plot checkboxes
        main_checks = (
            f.plot_avg_p,
            f.plot_base_p,
            f.plot_breech_p,
            f.plot_stag_p,
            f.plot_stag_l,
            f.plot_vel,
            f.plot_nozzle_v,
            f.plot_burnup,
            f.plot_eta,
        )
        for check in main_checks:
            check.trace_add("write", lambda *_: self.update_main_plot())

        # Aux plot checkboxes
        for check in (f.trace_hull, f.trace_press):
            check.trace_add("write", lambda *_: self.update_aux_plot())

        # Guide graph checkboxes
        for entry in (f.guide_plot_travel, f.guide_plot_volume, f.guide_plot_burnout, f.guide_chamber_ruler):
            entry.trace_add("write", lambda *_: self.update_guide_graph())

    def on_state_change(self, type_option: GunType) -> None:
        """Update plot checkbox visibility based on gun type."""
        f = self.frame

        if type_option == GunType.CONVENTIONAL:
            f.plot_nozzle_v.remove()
            f.plot_breech_p.localize("plotBreechP")
            f.plot_stag_p.remove()
            f.plot_stag_l.remove()
            f.plot_eta.remove()

        elif type_option == GunType.RECOILLESS:
            f.plot_nozzle_v.restore()
            f.plot_breech_p.localize("plotNozzleP")
            f.plot_stag_p.restore()
            f.plot_stag_p.localize("plotStagP")
            f.plot_stag_l.restore()
            f.plot_stag_l.localize("plotStagL")
            f.plot_eta.restore()
            f.plot_eta.localize("plotEtaEsc")

        else:
            raise ValueError

    @property
    def context(self):
        """Build matplotlib context from current theme colors."""
        style = ttk.Style(self.frame)
        bgc = str(style.lookup("TEntry", "background"))
        fgc = str(style.lookup("TEntry", "foreground"))
        fbgc = str(style.lookup("TEntry", "fieldbackground")) or bgc

        return {
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
            "xaxis.labellocation": "right",
            "yaxis.labellocation": "top",
            "figure.facecolor": bgc,
            "figure.edgecolor": fgc,
            "axes.edgecolor": fgc,
            "axes.labelcolor": fgc,
            "axes.facecolor": fbgc,
            "text.color": fgc,
            "xtick.color": fgc,
            "ytick.color": fgc,
        }

    @property
    def master(self) -> InteriorBallisticsFrame:
        """The top-level frame owning the simulation state."""
        return self.frame.master

    @property
    def theme_cmap(self):
        return mpl.colormaps[THEMES[Theme(self.frame.theme_name_var.get())]["cmap"]]

    def _setup_canvas(self, canvas):
        """Clear and prepare a canvas for plotting."""
        canvas.figure.clear()
        canvas.figure.set_facecolor(self.context["figure.facecolor"])
        return canvas.figure.add_subplot(111)

    def update_main_plot(self):
        """Main pressure/velocity plot."""
        with plt.rc_context(self.context):
            canvas = self.frame.fig_canvas
            ax = self._setup_canvas(canvas)
            ax_p, ax_v = ax.twinx(), ax.twinx()
            config, gun, gun_result = self.master.config, self.master.gun, self.master.gun_result
            if config and gun and gun_result:
                v_tgt = config.design_velocity
                p_tgt = config.design_pressure / 1e6  # Pa -> MPa to match axis
                gun_type = config.gun_type
                dom = config.domain

                xs, vs, vxs, pas, pss, pbs, p0s, psis, etas, stags = [], [], [], [], [], [], [], [], [], []

                for entry in gun_result.table_data:
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

                    if tag == config.pressure_control_point:
                        x_peak = (time * 1e3) if dom == Domain.TIME else travel
                        # noinspection PyTypeChecker
                        ax_p.spines.right.set_position(("data", x_peak))

                    if dom == Domain.TIME:
                        xs.append(time * 1000)
                    elif dom == Domain.LEN:
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

                if self.frame.plot_breech_p.get():
                    ax_p.plot(
                        xs,
                        pbs,
                        c="xkcd:goldenrod",
                        label=(
                            self.frame.get_loc_str("figBreech" if gun_type == GunType.CONVENTIONAL else "figNozzleP")
                        ),
                    )

                if gun_type == GunType.RECOILLESS:
                    if self.frame.plot_stag_p.get():
                        ax_p.plot(xs, p0s, "seagreen", label=self.frame.get_loc_str("figStagnation"))

                    if self.frame.plot_stag_l.get():
                        ax.plot(xs, stags, "mediumorchid", label=self.frame.get_loc_str("figStagL"))

                    if self.frame.plot_nozzle_v.get():
                        ax_v.plot(xs, vxs, "steelblue", label=self.frame.get_loc_str("figNozzleV"))

                    if self.frame.plot_eta.get():
                        ax.plot(xs, etas, "maroon", label=self.frame.get_loc_str("figOutflow"))

                if self.frame.plot_avg_p.get():
                    ax_p.plot(xs, pas, "tab:green", label=self.frame.get_loc_str("figAvgP"))

                if self.frame.plot_base_p.get():
                    ax_p.plot(xs, pss, "yellowgreen", label=self.frame.get_loc_str("figShotBase"))

                if gun_type in (GunType.CONVENTIONAL, GunType.RECOILLESS):
                    ax_p.axhline(float(p_tgt), c="tab:green", linestyle=":", label=self.frame.get_loc_str("figTgtP"))

                if self.frame.plot_vel.get():
                    ax_v.plot(xs, vs, "tab:blue", label=self.frame.get_loc_str("figShotVel"))
                ax_v.axhline(v_tgt, c="tab:blue", linestyle=":", label=self.frame.get_loc_str("figTgtV"))

                if self.frame.plot_burnup.get():
                    ax.plot(xs, psis, c="crimson", label=self.frame.get_loc_str("figPsi"))

                lines_labeled = []
                for lines, xvals in zip(
                    (ax_p.get_lines(), ax.get_lines(), ax_v.get_lines()),
                    (
                        (0.2 * xs[-1] + 0.8 * x_peak, xs[-1]),
                        (0, xs[-1]),
                        (x_peak, 0.2 * xs[-1] + 0.8 * x_peak),
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

                if dom == Domain.TIME:
                    ax.set_xlabel(self.frame.get_loc_str("figTimeDomain"))
                elif dom == Domain.LEN:
                    ax.set_xlabel(self.frame.get_loc_str("figLenDomain"))

            canvas.draw_idle()

    def update_aux_plot(self):
        """Auxiliary pressure trace plot."""
        with plt.rc_context(self.context):
            canvas = self.frame.aux_canvas
            aux_ax = self._setup_canvas(canvas)
            aux_ax_h = aux_ax.twinx()
            aux_ax_h.yaxis.tick_right()
            gun_result = self.master.gun_result
            if gun_result:
                p_trace = gun_result.pressure_trace
                cmap = self.theme_cmap
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

                    if self.frame.trace_press.get():
                        aux_ax.plot(x, y, c=color, alpha=alpha, ls=linestyle)

                aux_ax.set_xlim(left=0, right=x_max)
                aux_ax.set_ylim(bottom=0, top=y_max * 1.15)

                aux_ax.tick_params(axis="y", colors="tab:green")
                aux_ax_h.tick_params(axis="y", colors="tab:blue")
                aux_ax.tick_params(axis="x")

                aux_ax.set_xlabel(self.frame.get_loc_str("figAuxDomain"))

                aux_ax.yaxis.label.set_color("tab:green")
                aux_ax.set_ylabel("MPa")

                aux_ax_h.yaxis.set_ticks_position("right")
                aux_ax_h.yaxis.set_label_position("right")

                aux_ax_h.yaxis.label.set_color("tab:blue")
                aux_ax_h.set_ylabel("mm")

                h_trace = gun_result.outline

                if h_trace is not None and self.frame.trace_hull.get():
                    x_hull = list(entry.x for entry in h_trace)
                    r_in = list(entry.r_in * 1e3 for entry in h_trace)
                    r_out = list(entry.r_ex * 1e3 for entry in h_trace)
                    r_pej = list(entry.r_pej * 1e3 for entry in h_trace)

                    aux_ax_h.fill_between(
                        x_hull, r_in, r_pej, alpha=0.5 if self.frame.trace_press.get() else 0.8, color="tab:orange"
                    )
                    aux_ax_h.fill_between(
                        x_hull, r_pej, r_out, alpha=0.5 if self.frame.trace_press.get() else 0.8, color="tab:blue"
                    )

                    aux_ax.set_xlim(left=min(x_hull))

                aux_ax_h.set_ylim(bottom=0)

            canvas.draw_idle()

    def update_geom_plot(self):
        """Geometry σ(Z) plot."""
        with plt.rc_context(self.context):
            canvas = self.frame.geom_canvas
            geom_ax = self._setup_canvas(canvas)
            n = 100
            prop = self.master.prop
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
                geom_ax.yaxis.set_ticks([i * 0.5 for i in range(math.ceil(max(ys) / 0.5) + 1)])

            canvas.draw_idle()

    def update_guide_graph(self):
        """Guide graph with contours."""

        style = ttk.Style(self.frame)
        fgc = str(style.lookup("TFrame", "foreground"))
        cmap = self.theme_cmap.reversed()

        with plt.rc_context(self.context):
            canvas = self.frame.guide_canvas
            guide_ax = self._setup_canvas(canvas)
            gun, guide_results, config = self.master.gun, self.master.guide_results, self.master.config

            if guide_results:
                load_densities = list(line.load_density for line in guide_results.lines)
                charge_masses = list(line.charge_mass for line in guide_results.lines)

                delta_max = max(load_densities)
                delta_min = min(load_densities)
                w_max = max(charge_masses)
                w_min = min(charge_masses)
                max_cv = w_max / delta_min
                min_cv = w_min / delta_max

                left_diagonal_chamber_volume = w_min / delta_min
                right_diagonal_chamber_volume = w_max / delta_max

                if gun and config:
                    guide_ax.scatter(
                        config.charge_mass / config.chamber_volume,
                        config.charge_mass,
                        c=fgc,
                        marker="x",
                        s=FONTSIZE * 4,
                    )

                if self.frame.guide_chamber_ruler.get():

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

                guide_ax.set_xlabel(self.frame.get_loc_str("guideLDDomain"))
                guide_ax.set_ylabel(self.frame.get_loc_str("guideCMDomain"))

                titles = []

                for index, show, scaling, levels_func, title_loc_str, linestyle, unit in zip(
                    ("length_gun", "volume", "burnout"),
                    (
                        self.frame.guide_plot_travel.get(),
                        self.frame.guide_plot_volume.get(),
                        self.frame.guide_plot_burnout.get(),
                    ),
                    (1, 1000, 1),
                    (get_adaptive_scale, get_adaptive_scale, get_01_scale),
                    ("guideTravelTitle", "guideBVTitle", "guideBurnoutTitle"),
                    ("--", "-.", ":"),
                    (" m", " L", ""),
                ):
                    if not show:
                        continue

                    titles.append(self.frame.get_loc_str(title_loc_str))
                    results = list(line.__getattribute__(index) * scaling for line in guide_results.lines)
                    result_levels = levels_func(results)

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

            canvas.draw_idle()

    def get_figure(self, save_type: FigureType):
        """Get the matplotlib figure for the specified type."""
        if save_type == FigureType.MAIN:
            return self.frame.fig_canvas.figure
        elif save_type == FigureType.AUX:
            return self.frame.aux_canvas.figure
        elif save_type == FigureType.GUIDE:
            return self.frame.guide_canvas.figure
        elif save_type == FigureType.GEOM:
            return self.frame.geom_canvas.figure
        return None

    def export_figure(self, save_type: FigureType, file_name: str):
        """Export a figure as PNG with proper matplotlib context."""
        fig = self.get_figure(save_type)
        if fig is None:
            return False

        with plt.rc_context(self.context):
            fig.savefig(file_name, transparent=True, dpi=300)

        return True
