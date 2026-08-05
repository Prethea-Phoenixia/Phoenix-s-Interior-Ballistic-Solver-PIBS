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
- **Ballistics core**: `pibs/ballistics/` — `gun.py`, `recoilless.py`, `constrained*.py`, `num/` (RK45/UMF integrators), `prop/` (propellant models)
- **Resources**: loaded via `resolve_path()` from `pibs/misc.py` — handles both dev mode and frozen PyInstaller builds.

## Localization

- File: `pibs/ui/localization.json`
- After editing keys, run: `.venv\Scripts\python.exe sort_localization.py` (sorts keys for consistency)

## UI Widget Conventions

- `pibs/localized_widget.py` defines all localized widgets and `RowBuilder`.
- Widgets do **not** self-grid in `__init__`; they are placed via a `place(row, col, ...)` method called by `RowBuilder`.
- `RowBuilder` is the single source of truth for row placement — don't mix direct `.grid()` calls with RowBuilder-managed widgets.
- For checkboxes used as LabelFrame headers: use `widget.as_labelwidget()` instead of `.place()`.

## Packaging

- PyInstaller via `generate_executable.py`; spec file generated in root.
- Resources bundled via `--add-data` flags; always use `resolve_path()` for resource access.
