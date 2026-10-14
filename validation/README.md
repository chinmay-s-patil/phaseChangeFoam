# Natural Convection & Phase Change Validation Ladder

This directory contains the self-contained validation suite for the `phaseChangeFoam` suite of solvers (`phaseChangeMultiRegionFoam` and `phaseChangeMultiRegionSimpleFoam`).

It implements a progressive two-rung validation hierarchy verifying both pure natural convection buoyancy and coupled melting phase change against classic literature benchmarks.

---

## Validation Hierarchy

```
                    [ Rung 2: Gallium Melting Benchmark ]
                    Brent, Voller, Reid (1988)
                    Coupled thermal phase change (enthalpy-porosity) +
                    Darcy drag damping + Natural convection (Ra ~ 6e5, Pr = 0.0216)
                                   ▲
                                   │  (Builds on pure natural convection)
                                   │
                    [ Rung 1: Natural Convection Cavity ]
                    de Vahl Davis (1983)
                    Pure buoyancy circulation (Ra = 10^5, Pr = 0.71)
                    Steady-state solver: phaseChangeMultiRegionSimpleFoam
```

---

## Benchmarks Overview

### Rung 1: de Vahl Davis (1983) Natural Convection Cavity
- **Directory**: [`deVahlDavis2D/`](deVahlDavis2D/)
- **Solver**: `phaseChangeMultiRegionSimpleFoam` (steady-state Boussinesq)
- **Physics**: Laminar natural convection in a differentially heated square cavity ($L = 0.1\text{ m}$) at $\text{Ra} = 10^5, \text{Pr} = 0.71$.
- **Benchmark Metric**: Mean Nusselt number along the hot wall ($\overline{\text{Nu}} = 4.519$).
- **Status**: **PASSED** ($\overline{\text{Nu}} = 4.5830$, relative difference **1.42%**, tolerance $\le 5\%$).
- **Runtime**: $\approx 2\text{ s}$ (400 SIMPLE iterations).

### Rung 2: Brent, Voller, Reid (1988) Gallium Melting Benchmark
- **Directory**: [`gallium2D/`](gallium2D/)
- **Solver**: `phaseChangeMultiRegionFoam` (transient PIMPLE with enthalpy-porosity)
- **Physics**: Melting of pure gallium in a rectangular enclosure ($8.89\text{ cm} \times 6.35\text{ cm}$) subject to a constant hot wall temperature ($T_h = 311.15\text{ K}$) and initial/cold wall temperature ($T_c = 301.45\text{ K}$), with melting point $T_m = 302.93\text{ K}$.
  - Stefan number: $\text{Ste} = c_p \Delta T / L = 0.039$
  - Prandtl number: $\text{Pr} = \nu / \alpha = 0.0216$
  - Rayleigh number: $\text{Ra} = g \beta \Delta T H^3 / (\nu \alpha) \approx 6 \times 10^5$
- **Benchmark Metrics**:
  1. Liquid volume fraction $f_L(t)$ at $t = 2\text{ min}$ (120 s) against numerical benchmark.
  2. Melt front contour position $x_{front}(y)$ at $t = 2\text{ min}$ digitized from Figure 9 of Brent et al. (1988).
- **Status**: **PASSED** (Front position accurate to within $1.2$ grid cell widths across the full cavity width).
- **Runtime**: $\approx 35\text{ s}$ for 120 s physical simulation ($1,200$ time steps).

---

## Reference Data & Literature Sources
- **`papers/brent1988.pdf`**: Direct copy of Brent, Voller, Reid (1988), *Numerical Heat Transfer*, 13(3), pp. 297–318.
- **`gallium2D/reference_data/`**: Digitized melt front coordinates ($t = 2, 6, 10, 17\text{ min}$) and liquid fraction curve extracted from Brent et al. (1988).

---

## Running the Validation Ladder

Run the complete validation ladder from the repository root:

```bash
# Run the complete ladder (Rung 1 & Rung 2):
python3 validation/runAll.py

# Run fast gate (Rung 1 de Vahl Davis only, ~2 s):
python3 validation/runAll.py --fast

# Run individual benchmarks:
python3 validation/deVahlDavis2D/runTest.py
python3 validation/gallium2D/runTest.py --endTime 120 --plot
```
