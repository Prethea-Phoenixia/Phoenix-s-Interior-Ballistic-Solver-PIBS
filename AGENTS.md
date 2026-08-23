# PIBS — Agent Instructions

## Run / Test

- **Run app**: `.venv\Scripts\python.exe run_pibs.py` — always use project .venv. App is a tkinter GUI; run directly and watch terminal for errors. User will close window when done.
- **No test suite** — verification is manual (launch app, run calculation, save/load, theme switch).

## Lint / Format

- **black**: `line-length = 120` — run via `.venv\Scripts\python.exe -m black <files>`
- **isort**: `profile = "black"` — run via `.venv\Scripts\python.exe -m isort <files>`
- Pre-commit hook (`pre-commit`) also runs `sort_localization.py` and `generate_executable.py`, but isort/PyInstaller are often not in PATH; black is the only reliable step.

## Architecture

- **Entry**: `run_pibs.py` → `pibs.interior_ballistics.main()` → `PIBS(Tk)` → `InteriorBallisticsFrame`
- **UI layer**: `pibs/interior_ballistics_frame.py`, `notebook_frame.py`, `info_frame.py`, `table_frame.py` (tkinter + matplotlib)
- **Config**: `pibs/config.py` — `SimulationConfig` dataclass (all fields required, no defaults except `logger`). UI gathers values into this; dispatch passes it to ballistics core.
- **Ballistics core**: `pibs/ballistics/` — `gun.py`, `recoilless.py`, `constrained*.py`, `config.py` (typed dataclasses with `__post_init__` validation: `GunGeometry`, `PropellantLoad`, `Solver`, etc.), `num/` (RK45/UMF integrators), `prop/` (propellant models)
- **Dispatch**: `pibs/dispatch.py` — `calculate()` and `guide()` run in separate processes; `sim_config_to_ballistics()` converts `SimulationConfig` to ballistics objects (triggers validation).
- **Resources**: loaded via `resolve_path()` from `pibs/misc.py` — handles both dev mode and frozen PyInstaller builds.

## Considered Design Decisions

**"God objects" in ballistics core are intentional, not technical debt.**

- `Gun`/`Recoilless` own the full simulation pipeline (ODE setup, integration, peak finding, sampling, pressure traces). Splitting into separate `PressureCalculator`, `BurnModel`, `MotionSolver` classes would add indirection over inherently coupled state (`z`, `l_bar`, `v_bar`, `p_bar`) without real benefit for this scope.
- `InteriorBallisticsFrame` (1400+ lines) is a tkinter God frame — framework limitation, not design oversight.
- Refactors should not break these classes apart "for purity." Tight coupling here tracks the physics, not accidental complexity.

**Input validation lives in the domain layer, not the UI.**

- Ballistics config dataclasses use `__post_init__` to validate constraints (positive values, ranges, membership).
- The UI performs only **normalization** (empty→0.0, str→float via focus-out formatters), not constraint checking.
- Invalid values surface as `ValueError` at `dispatch.sim_config_to_ballistics()` — caught by the existing error handling infrastructure.
- This ensures consistent validation across all entry points (GUI, CLI, API, tests).
- Valid value collections (`VALID_SOLUTION_METHODS`, `VALID_PRESSURE_POINTS`, etc.) are defined in `ballistics/__init__.py` alongside their `Literal` type aliases.

## Save / Load

- Save files use **widget descriptive strings** as JSON keys (from `loc.get_descriptive()`), NOT magic keys.
- Load matches JSON keys to widgets via `loc_dict`; no round-trip through `SimulationConfig`.

## Localization

- File: `pibs/ui/localization.json`
- After editing keys, run: `.venv\Scripts\python.exe sort_localization.py` (sorts keys for consistency)

## UI Widget Conventions

- `pibs/localized_widget.py` defines all localized widgets and `RowBuilder`.
- Widgets do **not** self-grid in `__init__`; they are placed via a `place(row, col, ...)` method called by `RowBuilder`.
- `RowBuilder` is the single source of truth for row placement — don't mix direct `.grid()` calls with RowBuilder-managed widgets.
- For checkboxes used as LabelFrame headers: use `widget.as_labelwidget()` instead of `.place()`.

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
- Before committing a refactor that touches ballistics code (`gun.py`, `recoilless.py`, `constrained*.py`), compare against the previous version to ensure no technical docs were lost.

## Packaging

- PyInstaller via `generate_executable.py`; spec file generated in root.
- Resources bundled via `--add-data` flags; always use `resolve_path()` for resource access.
