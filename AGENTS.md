# PIBS — Agent Instructions

## Run / Test

- **Run app**: `.venv\Scripts\python.exe run_pibs.py` — always use project .venv. App is a tkinter GUI; run directly and
  watch terminal for errors. User will close window when done.
- **Unit tests**: `python -m unittest discover -s tests` (numerical primitives only; set PYTHONPATH to repo root).
  Everything else is verified manually (launch app, run calculation, save/load, theme switch).
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

## Architecture

```
run_pibs.py
    │
    ▼
PIBS (tk.Tk, pibs.py)
    ├─ InteriorBallisticsFrame (ib_frame.py)   input widgets; .config property builds SimConfig
    ├─ NotebookFrame (nb_frame.py)
    │    ├─ InfoFrame    (info_frame.py)    reads SimResult.params
    │    ├─ TableFrame   (tbl_frame.py)     reads SimResult.table
    │    └─ PlotManager  (plot_manager.py)  reads SimResult.table / .outline
    └─ file_io.py · mode_manager.py · localized.py · theme.py · tip.py
    │
    │  SimConfig pickled through multiprocessing.Queue
    ▼
run.calculate(cfg)  (spawned worker; QueueHandler logging)
    ├─ if cfg.constrained: solve_constrained() / find_min_v()   (constrained.py)
    │      writes solved web_thickness / barrel_length / chamber_volume back into cfg
    ├─ simulate(cfg) -> SimResult                               (simulate.py)
    └─ if cfg.compute_guide: guide_graph(cfg) -> GuideResults   (guidegraph.py, Pool)

pibs/ballistics
│
│  data structs (no behavior)            functions (the pipeline)
│  ──────────────────────────            ──────────────────────────────────
│  SimConfig, Material  (config.py) ──►  derive_params(cfg) -> Params      (params.py)
│  Params               (params.py, frozen)        │
│  SimResult, TableEntry, PressureTraceEntry,      ▼
│    PressureProbePoint, OutlineEntry          simulate(cfg) -> SimResult  (simulate.py)
│    (generics.py)                               ├─ _gun_model / _recoilless_model
│  GuideResults           (guidegraph.py)        │   (pressure closures + row building)
│  enums (__init__.py): Domain, Point, GunType,  ├─ _integrate (rkf driver)
│    SolutionMethod, OptimizationTarget          └─ barrel_monoblock / barrel_autofrettage,
│                                                │   gun_structure / recoilless_structure
│  odes.py                                       │
│    make_gun_rhs / make_recoilless_rhs          nozzle.py
│      master time-domain RHS (canonical state     get_cf / get_ar / get_pr
│      (t,z,l,v,η,τ); guns carry dummy η=0,τ=1)    (nozzle theory, Hunt 1953)
│    in_domain(rhs, T|Z|L|V) — chain-rule
│      re-parametrization; shared by simulate
│      and constrained
│
│  constrained.py
│    solve_constrained() ─┬─ _solve_gun          (iterates chambrage correction
│                         └─ _solve_recoilless    in a loop; reuses odes.py rhs)
│    find_min_v() · maximum_load_fraction()
│    validate_load_fraction() · get_f()
│    minimum_load_fraction()
│
│  prop/ (real polymorphism)                 num/ (numerical primitives)
│    Geometry (ABC)                            rkf     adaptive RK integrator
│      ├─ SimpleGeometry (Enum)                dekker  bracketed root finder,
│      │    SPHERE CYLINDER TUBE STRIP         │      single relative tol
│      └─ MultPerfGeometry (Enum)              gss     golden-section search
│           SEVEN_PERF_CYLINDER, *_ROSETTE,    utils   record merge, find_last_le
│           NINETEEN_PERF_*, ...
│           SEVEN_PERF_CYLINDER, *_ROSETTE, NINETEEN_PERF_*, ...
│    Composition   propellant thermochemistry (f, α, ρ_p, θ, u₁, n; propellants.csv)
│    Propellant    Composition + Geometry + grain dims -> form functions ψ(z), σ(z);
│                  mixed charges weight multiple geometries
```

## Design Decisions

**The ballistics core is functional, not object-oriented.**

- One flat config: `SimConfig` (`pibs/ballistics/config.py`) — all simulation inputs in SI, validated in
  `__post_init__`. Do NOT reintroduce per-concern config classes (GunConfig/NozzleConfig/...) or a UI-side
  aggregate with converter properties; each parameter has exactly one name and one home.
- `params.derive_params(cfg) -> Params` computes all ~20 coupled derived scalars once; the frozen `Params`
  dataclass is the shared namespace the ODE closures and the UI (`info_frame`, `plot_manager`) read from.
  Keep it a plain data struct — no behavior.
- `simulate.simulate(cfg) -> SimResult` is the single pipeline for conventional and recoilless guns
  (conventional states carry dummy η=0, τ=1 so both share one code path). The ODE right-hand sides are
  defined once per gun type as time-domain masters in `odes.py`; integration domains (t/z/l/v) are
  chain-rule re-parametrizations via `in_domain()` — do not hand-write per-domain ODE variants.
  Constrained design lives in `constrained.py` as plain functions (`solve_constrained`, `find_min_v`,
  `validate_load_fraction`) reusing the same master RHS factories.
- `Propellant`/`Composition`/`Geometry` keep real polymorphism (mixed-charge form functions, geometry
  dispatch) — that's where classes earn their keep. `InteriorBallisticsFrame` (~1000 lines) remains a
  tkinter God frame — framework limitation, not design oversight.

**Input validation lives in the domain layer, not the UI.**

- `SimConfig.__post_init__` validates constraints (design fields only when `constrained`, structural only
  when a material is set).
- The UI performs only **normalization** (empty→0.0, str→float via focus-out formatters), not constraint checking.
- Invalid values surface as `ValueError` from `SimConfig` construction or `derive_params` — caught by existing
  error handling (`run.calculate` logs and returns None results).

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
- Before committing a refactor that touches ballistics code (`params.py`, `simulate.py`, `nozzle.py`,
  `constrained.py`), compare against the previous version to ensure no technical docs were lost.

## Gotchas

**Logging**: Python's `logging.lastResort` handler writes unformatted messages to stderr when no handlers are found.
In pool children, setting the `pibs` logger level to CRITICAL is required — merely clearing handlers is insufficient.

**PlotManager**: Each plot update method must wrap plotting code in `plt.rc_context(self.context)` to apply theme/font
settings. `draw_idle()` must be called inside this context.

**Resources**: Always use `resolve_path()` from `pibs/misc.py` for resource access — handles both dev mode and frozen
PyInstaller builds.

## Packaging

- PyInstaller via `generate_executable.py`; spec file generated in root.
