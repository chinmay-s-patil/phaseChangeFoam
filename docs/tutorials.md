# Tutorials & Case Setup

> [!NOTE]
> **Work in Progress**
> 
> Standalone tutorial cases are currently not included in this repository. Tutorial cases and demonstration workflows are under active development and will be provided in a future release.

---

## Verifying the Solvers

In the meantime, verified simulation cases and automated runners are available in the test suites:

* **Verification Suite**: [`verification/`](../verification/) — Quantitative verification cases covering phase change kinetics, energy conservation, and thermal hysteresis.
* **Validation Benchmarks**: [`validation/`](../validation/) — Physical validation cases including pure gallium melting benchmarks.

For theoretical formulations, model equations, and keyword dictionaries, please refer to:
* [`docs/models.md`](models.md) — Governing equations, thermodynamic consistency, species evaporation, and PCM formulations.
* [`docs/dictionary.md`](dictionary.md) — Keyword reference tables and solver control templates.
* [`docs/validation.md`](validation.md) — Benchmark comparison procedures and quantitative convergence data.
