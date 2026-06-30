from __future__ import annotations

from .ballistics import POINT_BURNOUT, POINT_EXIT, POINT_PEAK_BREECH, POINT_PEAK_AVG, POINT_PEAK_SHOT, Recoilless, Gun
from .ballistics.gun import GunResult
from .ballistics.recoilless import RecoillessResult
from .localized_widget import LocalizedFrame
from .misc import to_si, round_sig, format_mass


class InfoFrame(LocalizedFrame):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.columnconfigure(0, weight=1)
        self.rowconfigure(0, weight=1)

        par_frm = self.add_localized_label_frame(self, label_loc_key="parFrmLabel")
        par_frm.grid(row=0, column=0, sticky="nsew")
        par_frm.columnconfigure(0, weight=1)

        i = 0
        self.ammo, i = (self.add_localized_2_line_display(parent=par_frm, row=i, label_loc_key="ammoLabel"), i + 2)
        self.pp, i = (
            self.add_localized_3_line_display(parent=par_frm, row=i, label_loc_key="ppLabel", tooltip_loc_key="ppText"),
            i + 3,
        )
        self.bop, i = self.add_localized_2_line_display(parent=par_frm, row=i, label_loc_key="bopLabel"), i + 2
        self.lx, i = (
            self.add_localized_3_line_display(
                parent=par_frm, row=i, label_loc_key="lxLabel", tooltip_loc_key="calLxText"
            ),
            i + 3,
        )
        self.mv, i = self.add_localized_2_line_display(parent=par_frm, row=i, label_loc_key="mvLabel"), i + 2
        self.va, i = (
            self.add_localized_2_line_display(
                parent=par_frm, row=i, label_loc_key="vaLabel", tooltip_loc_key="vinfText"
            ),
            i + 2,
        )
        self.te, i = (
            self.add_localized_2_line_display(
                parent=par_frm, row=i, label_loc_key="teffLabel", tooltip_loc_key="teffText"
            ),
            i + 2,
        )
        self.be, i = (
            self.add_localized_2_line_display(
                parent=par_frm, row=i, label_loc_key="beffLabel", tooltip_loc_key="beffText"
            ),
            i + 2,
        )
        self.pe, i = (
            self.add_localized_2_line_display(
                parent=par_frm, row=i, label_loc_key="peffLabel", tooltip_loc_key="peffText"
            ),
            i + 2,
        )
        self.pa, i = self.add_localized_2_line_display(parent=par_frm, row=i, label_loc_key="paLabel"), i + 2
        self.gm, i = self.add_localized_2_line_display(parent=par_frm, row=i, label_loc_key="gmLabel"), i + 2
        self.sj, i = self.add_localized_2_line_display(parent=par_frm, row=i, label_loc_key="sjLabel"), i + 2

        self.ld, i = (self.add_localized_2_line_display(parent=par_frm, row=i, label_loc_key="ldLabel"), i + 2)

        self.lf, i = (self.add_localized_2_line_display(parent=par_frm, row=i, label_loc_key="ldfLabel"), i + 2)

        par_frm.rowconfigure(i + 1, weight=1)

    def update_stats(self, gun: Gun | Recoilless, gun_result: GunResult | RecoillessResult, acc_exp: int):
        for entry in (
            *(self.te, self.be, self.pe, self.va, self.lx, self.ammo, self.pa, self.gm, self.pp, self.mv, self.bop),
            *(self.sj, self.ld, self.lf),
        ):
            entry.reset()

        caliber = gun.caliber
        eta_t, eta_b, eta_p = gun_result.get_eff()
        self.te.set(f"{eta_t * 100:.2f} %")
        self.be.set(f"{eta_b * 100:.2f} %")
        self.pe.set(f"{eta_p * 100:.2f} %")
        self.va.set(to_si(gun.v_j, unit="m/s"))
        self.lx.set(
            (
                f"{gun.l_g / caliber:.0f}" + " " + self.get_loc_str("calLabel"),
                f"{(gun.l_g + gun.l_c) / caliber:.0f}" + " " + self.get_loc_str("calLabel"),
            )
        )
        self.ammo.set(to_si(gun.l_c, unit="m"))
        ps = gun_result.read_table_data(POINT_PEAK_SHOT).shot_pressure
        self.pa.set(to_si(ps * gun.s / gun.m, unit="m/s²"))
        if gun_result.tube_mass:
            self.gm.set(format_mass(gun_result.tube_mass))
            self.gm.restore()
        else:
            self.gm.remove()

        peak_average_entry = gun_result.read_table_data(POINT_PEAK_AVG)
        peak_breech_entry = gun_result.read_table_data(POINT_PEAK_BREECH)
        self.pp.set(
            (
                f"{to_si(peak_average_entry.avg_pressure, unit='Pa'):}" + self.get_loc_str("mean"),
                f"{to_si(peak_breech_entry.breech_pressure, unit='Pa'):}" + self.get_loc_str("breech"),
            )
        )
        muzzle_entry = gun_result.read_table_data(POINT_EXIT)
        self.mv.set(to_si(muzzle_entry.velocity, unit="m/s"))
        try:
            burnout_entry = gun_result.read_table_data(POINT_BURNOUT)
            self.bop.set(f"{burnout_entry.travel / gun.l_g * 1e2:.2f} %")
        except ValueError:
            self.bop.set(self.get_loc_str("uncontained"))

        if isinstance(gun, Recoilless):
            self.sj.set(f"{to_si(gun.s_j, unit='m²', unit_dim=2):}")
            self.sj.restore()
        else:
            self.sj.remove()

        sigfig = acc_exp + 1
        w = gun.w
        cv = gun.vol_0

        rho = gun.rho_p
        self.lf.set(f"{round_sig(w / cv / rho * 100, n=sigfig)} %")
        self.ld.set(f"{round_sig(w / cv, n=sigfig)} kg/m³")
