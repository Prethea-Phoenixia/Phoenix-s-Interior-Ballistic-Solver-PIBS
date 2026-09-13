from __future__ import annotations

from dataclasses import dataclass

from .ballistics import (
    CONVENTIONAL,
    POINT_PEAK_AVG,
    POINT_PEAK_BREECH,
    POINT_PEAK_SHOT,
    POINT_PEAK_STAG,
    RECOILLESS,
    SOL_LAGRANGE,
)

MODE_FREE = "free"
MODE_CONSTRAINED = "constrained"
MODE_LOCK_LG = "lock_lg"
MODE_OPT = "opt"


@dataclass
class ModeState:
    """Represents the current UI mode and its constraints."""

    gun_type: str
    mode: str
    is_conventional: bool
    is_recoilless: bool


class ModeManager:
    """Manages widget state based on gun type and simulation mode."""

    def __init__(self, frame):
        self.frame = frame

    def get_mode_state(self) -> ModeState:
        """Determine current mode from UI state."""
        f = self.frame
        gun_type = f.get_gun_type()

        mode = MODE_FREE
        if f.is_constrained():
            mode = MODE_CONSTRAINED
        if f.is_optimization():
            mode = MODE_OPT
        if f.is_lock_length():
            mode = MODE_LOCK_LG

        return ModeState(
            gun_type=gun_type,
            mode=mode,
            is_conventional=gun_type == CONVENTIONAL,
            is_recoilless=gun_type == RECOILLESS,
        )

    def apply_gun_type(self, state: ModeState):
        """Apply gun-type-specific widget states."""
        f = self.frame

        if state.is_conventional:
            f.drop_gradient.enable()
            f.nozz_exp.disable()
            f.nozz_eff.disable()
            f.p_control.reset(
                {p: p for p in (POINT_PEAK_AVG, POINT_PEAK_SHOT, POINT_PEAK_BREECH)},
                overwrite=False,
            )

        elif state.is_recoilless:
            f.drop_gradient.set_by_obj(SOL_LAGRANGE)
            f.drop_gradient.disable()
            f.nozz_exp.enable()
            f.nozz_eff.enable()
            f.p_control.reset(
                {p: p for p in (POINT_PEAK_AVG, POINT_PEAK_SHOT, POINT_PEAK_STAG, POINT_PEAK_BREECH)},
                overwrite=False,
            )

    def apply_mode(self, state: ModeState):
        """Apply mode-specific widget states."""
        f = self.frame

        states = {
            f.v_tgt: False,
            f.p_tgt: False,
            f.opt: False,
            f.lock_Lg: False,
            f.min_web: False,
            f.lg_max: False,
            f.p_control: False,
            f.drop_opt_tgt: False,
            f.tbl_mm: True,
            f.web_mm: True,
        }

        if state.mode == MODE_CONSTRAINED:
            states.update(
                {
                    f.v_tgt: True,
                    f.p_tgt: True,
                    f.opt: True,
                    f.lock_Lg: True,
                    f.min_web: True,
                    f.lg_max: True,
                    f.p_control: True,
                    f.tbl_mm: False,
                    f.web_mm: False,
                }
            )
        elif state.mode == MODE_LOCK_LG:
            states.update(
                {
                    f.p_tgt: True,
                    f.lock_Lg: True,
                    f.min_web: True,
                    f.lg_max: True,
                    f.p_control: True,
                    f.tbl_mm: True,
                    f.web_mm: False,
                }
            )
        elif state.mode == MODE_OPT:
            states.update(
                {
                    f.v_tgt: True,
                    f.p_tgt: True,
                    f.opt: True,
                    f.min_web: True,
                    f.lg_max: True,
                    f.p_control: True,
                    f.drop_opt_tgt: True,
                    f.tbl_mm: False,
                    f.web_mm: False,
                }
            )

        for widget, enabled in states.items():
            if enabled:
                widget.enable()
            else:
                widget.disable()

        if state.mode in (MODE_CONSTRAINED, MODE_OPT) and state.is_conventional:
            f.max_iter.enable()
        else:
            f.max_iter.disable()

    def apply_groups(self):
        """Apply checkbox-group-based enabling."""
        f = self.frame
        self._enable_group(
            f.use_aux_grain,
            f.aux_grain_r1,
            f.aux_grain_r2,
            f.aux_web_ratio,
            f.aux_mass_ratio,
            f.aux_geom,
        )
        self._enable_group(
            f.use_material,
            f.material_yield,
            f.material_density,
            f.material_ssf,
            f.material_is_af,
        )
        self._enable_group(f.use_combustible, f.combustible_mass_kg, f.combustible_force_kJ__kg)

    @staticmethod
    def _enable_group(control, *widgets):
        for w in widgets:
            w.enable() if control.get() else w.disable()

    def update(self):
        """Main update method - call this on any state change."""
        state = self.get_mode_state()
        self.frame.notebook_frame.on_state_change(type_option=state.gun_type)
        self.apply_gun_type(state)
        self.apply_mode(state)
        self.apply_groups()
