from __future__ import annotations

from .ballistics import POINT_BURNOUT, POINT_EXIT, POINT_PEAK_AVG, POINT_PEAK_BREECH, POINT_PEAK_SHOT
from .ballistics.gun import Gun, GunResult
from .ballistics.recoilless import Recoilless, RecoillessResult
from .localized import LocalizedFrame, RowBuilder
from .misc import format_mass, round_sig, to_si


class InfoFrame(LocalizedFrame):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.columnconfigure(0, weight=1)
        self.rowconfigure(0, weight=1)

        par_frm = self.add_localized_label_frame(self, label_loc_key="parFrmLabel")
        par_frm.grid(row=0, column=0, sticky="nsew")
        par_frm.columnconfigure(0, weight=1)

        b = RowBuilder(self, par_frm)
        self.ammo = b.display_2(label_loc_key="ammoLabel")
        self.pp = b.display_3(label_loc_key="ppLabel", tooltip_loc_key="ppText")
        self.bop = b.display_2(label_loc_key="bopLabel")
        self.lx = b.display_3(label_loc_key="lxLabel", tooltip_loc_key="calLxText")
        self.mv = b.display_2(label_loc_key="mvLabel")
        self.va = b.display_2(label_loc_key="vaLabel", tooltip_loc_key="vinfText")
        self.te = b.display_2(label_loc_key="teffLabel", tooltip_loc_key="teffText")
        self.be = b.display_2(label_loc_key="beffLabel", tooltip_loc_key="beffText")
        self.pe = b.display_2(label_loc_key="peffLabel", tooltip_loc_key="peffText")
        self.pa = b.display_2(label_loc_key="paLabel")
        self.gm = b.display_2(label_loc_key="gmLabel")
        self.sj = b.display_2(label_loc_key="sjLabel")
        self.ld = b.display_2(label_loc_key="ldLabel")
        self.lf = b.display_2(label_loc_key="ldfLabel")

        par_frm.rowconfigure(b.current_row + 1, weight=1)

    def update_stats(self, gun: Gun | Recoilless | None, gun_result: GunResult | RecoillessResult, acc_exp: int):
        for entry in (
            *(self.te, self.be, self.pe, self.va, self.lx, self.ammo, self.pa, self.gm, self.pp, self.mv, self.bop),
            *(self.sj, self.ld, self.lf),
        ):
            entry.reset()

        if not gun:
            return

        caliber = gun.caliber
        eta_t, eta_b, eta_p = gun_result.get_eff()
        self.te.set(f"{eta_t * 100:.2f} %")
        self.be.set(f"{eta_b * 100:.2f} %")
        self.pe.set(f"{eta_p * 100:.2f} %")
        self.va.set(to_si(gun.v_j, unit="m/s"))
        self.lx.set(
            (
                f"{gun.l_g / caliber:.0f}" + " " + self.get_loc_str("calLabel"),
                f"{float((gun.l_g + gun.l_c) / caliber):.0f}" + " " + self.get_loc_str("calLabel"),
            )
        )
        self.ammo.set(to_si(gun.l_c, unit="m"))
        ps = gun_result.read_table_data(POINT_PEAK_SHOT).shot_pressure
        self.pa.set(to_si(ps * gun.s / gun.m, unit="m/s²"))

        self.gm.set(format_mass(gun_result.tube_mass) if gun_result.tube_mass else "N/A")

        peak_average_entry = gun_result.read_table_data(POINT_PEAK_AVG)
        peak_breech_entry = gun_result.read_table_data(POINT_PEAK_BREECH)
        self.pp.set(
            (
                f"{to_si(peak_average_entry.avg_pressure, unit='Pa')}" + self.get_loc_str("mean"),
                f"{to_si(peak_breech_entry.breech_pressure, unit='Pa')}" + self.get_loc_str("breech"),
            )
        )
        muzzle_entry = gun_result.read_table_data(POINT_EXIT)
        self.mv.set(to_si(muzzle_entry.velocity, unit="m/s"))
        try:
            burnout_entry = gun_result.read_table_data(POINT_BURNOUT)
            self.bop.set(f"{burnout_entry.travel / gun.l_g * 1e2:.2f} %")
        except ValueError:
            self.bop.set(self.get_loc_str("uncontained"))

        self.sj.set(f"{to_si(gun.s_j, unit='m²', unit_dim=2)}" if isinstance(gun, Recoilless) else "N/A")

        sigfig = acc_exp + 1
        w = gun.w
        cv = gun.vol_0

        rho = gun.rho_p
        self.lf.set(f"{round_sig(w / cv / rho * 100, n=sigfig)} %")
        self.ld.set(f"{round_sig(w / cv, n=sigfig)} kg/m³")
