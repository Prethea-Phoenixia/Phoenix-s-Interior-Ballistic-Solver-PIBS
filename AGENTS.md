# PIBS — Agent Instructions

## Run / Test

- **Run app**: `.venv\Scripts\python.exe run_pibs.py` — always use project .venv. App is a tkinter GUI; run directly and
  watch terminal for errors. User will close window when done.
- **No test suite** — verification is manual (launch app, run calculation, save/load, theme switch).
- **In-process GUI testing (welcome method)** — construct the real window without `mainloop()` to verify behavior
  headlessly:
  ```python
  root = PIBS(loc="English", debug=False)
  root.update()                          # event-loop pass; required before map/grid state is valid
  f.main_geom.set_by_obj(geom); f.update_geom(); root.update()   # same path a user click takes
  assert not f.grain_r1.input_widget.winfo_ismapped()             # assert real UI state, not just "no crash"
  root.destroy()
  ```
  Useful assertions: `winfo_ismapped()` (visibility), `cget("text")` (labels), `menu.entrycget(i, "label")` (menus),
  `frame.get_save_data()` (save keys). Run with `-X utf8` and wrap `sys.stdout` in a UTF-8 `TextIOWrapper` (console is
  GBK). Covers behavior, not visual layout — the user still eyeballs layout changes in the live window.

## Lint / Format

- **black**: `line-length = 120` — run via `.venv\Scripts\python.exe -m black <files>`
- **isort**: `profile = "black"` — run via `.venv\Scripts\python.exe -m isort <files>`
- Pre-commit hook (`pre-commit`) runs `black .`, `isort .`, and `sort_localization.py`. isort is often not in PATH here;
  black is the most reliable step.

## Design Decisions

**"God objects" in ballistics core are intentional, not technical debt.**

- `Gun`/`Recoilless` own the full simulation pipeline (ODE setup, integration, peak finding, sampling, pressure traces).
  Splitting into separate classes would add indirection over inherently coupled state without real benefit.
- `InteriorBallisticsFrame` (~1000 lines) is a tkinter God frame — framework limitation, not design oversight.
- Refactors should not break these classes apart "for purity." Tight coupling here tracks the physics.

**Input validation lives in the domain layer, not the UI.**

- Ballistics config dataclasses use `__post_init__` to validate constraints.
- The UI performs only **normalization** (empty→0.0, str→float via focus-out formatters), not constraint checking.
- Invalid values surface as `ValueError` at `dispatch.sim_config_to_ballistics()` — caught by existing error handling.

## Conventions

**Save / Load**

- Save files use **widget descriptive strings** as JSON keys (from `loc.get_descriptive()`), NOT magic keys.

**Localization**

- After editing keys in `pibs/ui/localization.json`, run: `.venv\Scripts\python.exe sort_localization.py`

**UI Widgets**

- Widgets **self-grid in `__init__`** using `parent`/`row`/`col` passed via kwargs. `RowBuilder` injects `parent` and the
  current `row`, then increments its counter — there is no separate `place()` method.
- `RowBuilder` is the single source of truth for row placement — don't hand-manage row indices for builder-managed widgets.
- For a checkbox used as a LabelFrame header: create it with `skip_grid=True`, then pass `labelwidget=check.check_widget`
  to `add_localized_label_frame`.

## Documentation Policy

**CRITICAL: Preserve technical/research documentation.**

- Code docstrings that describe API usage can be updated or removed during refactors.
- **Technical documentation must never be deleted** — this includes:
    - Textbook/paper references (e.g., 金 2014, Hunt 1953, 鲍廷钰 1995)
    - Mathematical formulas and derivations
    - Physics explanations (ODE domains, pressure conversions, nozzle theory)
    - ASCII diagrams explaining geometry or flow
    - Citations to specific equations or page numbers
- When refactoring, **move** technical docs to their new location — never drop them.
- Before committing a refactor that touches ballistics code (`gun.py`, `recoilless.py`, `constrained*.py`), compare
  against the previous version to ensure no technical docs were lost.

## Gotchas

**Logging**: Python's `logging.lastResort` handler writes unformatted messages to stderr when no handlers are found.
In pool children, setting the `pibs` logger level to CRITICAL is required — merely clearing handlers is insufficient.

**PlotManager**: Each plot update method must wrap plotting code in `plt.rc_context(self.context)` to apply theme/font
settings. `draw_idle()` must be called inside this context.

**Resources**: Always use `resolve_path()` from `pibs/misc.py` for resource access — handles both dev mode and frozen
PyInstaller builds.

## Packaging

- PyInstaller via `generate_executable.py`; spec file generated in root.
