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
  1. Liquid volume fraction evolution $f_L(t)$ evaluated across the full $17\text{ min}$ ($1020\text{ s}$) benchmark ($t = 2, 6, 10, 17\text{ min}$).
  2. Melt front contour position $x_{front}(y)$ at $t = 2, 6, 10, 17\text{ min}$ compared against numerical and experimental front positions.
- **Status**: **UNDER RIGOROUS VALIDATION / AUDIT** (Comprehensive spatial mesh, time-step convergence, and exact per-step energy balance audit).
- **Runtime**: $\approx 35\text{ s}$ for 120 s short gate; $\approx 4\text{--}5\text{ min}$ for full 17 min simulation ($10,200$ time steps at $\Delta t = 0.1\text{ s}$).

---

## Reference Data & Literature Sources
- **Brent, Voller, Reid (1988)**: "Enthalpy-Porosity Technique for Modeling Convection-Diffusion Phase Change: Application to the Melting of a Pure Metal", *Numerical Heat Transfer*, 13(3), pp. 297–318. DOI: [10.1080/10407788808913615](https://doi.org/10.1080/10407788808913615).
- **Hannoun, Alexiades, Mai (2003)**: "Resolving the Controversy over Tin and Gallium Melting in a Rectangular Cavity", *Numerical Heat Transfer, Part B: Fundamentals*, 44(3), pp. 253–276. DOI: [10.1080/716100520](https://doi.org/10.1080/716100520).
- **Ye & Arıcı (2023)**: "Redefined-interface-error-based 2D verification and validation study on natural convection melting: Benchmark with pure Gallium", *International Communications in Heat and Mass Transfer*, 147, 106952. DOI: [10.1016/j.icheatmasstransfer.2023.106952](https://doi.org/10.1016/j.icheatmasstransfer.2023.106952).
- **Gau & Viskanta (1986)**: "Melting and Solidification of a Metal System—I. Melting from a Vertical Wall", *Int. J. Heat Mass Transfer*, 29(12), pp. 1951–1963.
- **`gallium2D/reference_data/`**: Digitized melt front coordinates ($t = 2, 6, 10, 17\text{ min}$) and liquid fraction curves from Brent et al. (1988) and subsequent literature.

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
