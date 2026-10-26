# phaseChangeFoam

[![CI](https://github.com/chinmay-s-patil/phaseChangeFoam/actions/workflows/ci.yml/badge.svg)](https://github.com/chinmay-s-patil/phaseChangeFoam/actions/workflows/ci.yml)
[![OpenFOAM](https://img.shields.io/badge/openfoam-v2412-informational)](https://www.openfoam.com)
[![Docker](https://github.com/chinmay-s-patil/phaseChangeFoam/actions/workflows/docker.yml/badge.svg)](https://github.com/chinmay-s-patil/phaseChangeFoam/actions/workflows/docker.yml)
[![Docker image](https://img.shields.io/badge/docker%20image-ghcr.io%2Fchinmay--s--patil%2Fphasechangefoam-0085a0?logo=docker&logoColor=white)](https://github.com/chinmay-s-patil/phaseChangeFoam/pkgs/container/phasechangefoam)
[![Release](https://img.shields.io/github/v/release/chinmay-s-patil/phaseChangeFoam?color=blue)](https://github.com/chinmay-s-patil/phaseChangeFoam/releases)
[![License](https://img.shields.io/badge/license-GPL--3.0-green.svg)](LICENSE)

`phaseChangeFoam` is a multi-region conjugate heat transfer (CHT) solver suite for OpenFOAM with integrated phase change physics. It provides solid-liquid phase change materials (PCM) with thermal hysteresis and Darcy mushy-zone drag, as well as species evaporation and boiling kinetics across coupled solid and fluid domains.

> **Status**: Research code tested on OpenFOAM-v2412. Release v0.1.0; interfaces and keyword specifications are subject to refinement.

---

## Solvers and Models

### Solvers
| Solver | Description | Time Formulation |
| :--- | :--- | :--- |
| `phaseChangeMultiRegionFoam` | Multi-region transient CHT with species phase change and solid-liquid PCM | PIMPLE |
| `phaseChangeMultiRegionSimpleFoam` | Multi-region steady-state CHT with natural convection and thermal phase change | SIMPLE |

### Phase Change Models
| Model | Domain | Formulation | Primary Mechanisms |
| :--- | :--- | :--- | :--- |
| `EHCDarcy` | Solid / Fluid | Effective Heat Capacity + Darcy drag | Matrix-implicit capacity, Carman-Kozeny velocity damping |
| `EHC` | Solid / Fluid | Effective Heat Capacity | Path-integral capacity, subcooling hysteresis (no velocity damping) |
| `enthalpyPorosity` | Solid / Fluid | Enthalpy-porosity (Voller and Prakash) | Newton projection on liquid fraction, mushy-zone Darcy drag |
| `Lee` | Fluid | Volumetric species mass/energy transfer | Evaporation and condensation, optional Clausius-Clapeyron Tsat(p) |
| `constantSource` | Fluid | Fixed-rate mass/energy source | Verification and testing baseline |

---

## Differences from chtMultiRegionFoam

`phaseChangeFoam` is derived from OpenFOAM's `chtMultiRegionFoam` (OpenCFD / ESI-OpenCFD), inheriting its multi-region conjugate heat transfer structure and monolithic coupled energy framework (`coupledEnergyField`, `useImplicit`). The primary additions are:

* **Thermal hysteresis in solid-liquid PCM**: Stateful direction-dependent melting and freezing paths with subcooling tracking (`phaseState`, `T_reversal`).
* **Effective Heat Capacity formulation (`EHCDarcy`, `EHC`)**: Directly incorporates latent heat into cell capacity $C_{p,\text{eff}}$, avoiding the explicit source term under-relaxation oscillations characteristic of enthalpy-porosity.
* **Adaptive time-step limiting (`maxAlphaCo`)**: Constrains $\Delta t$ based on maximum cell-wise rate of phase fraction change, preventing the interface from stepping over cells without resolving the transition.
* **Volumetric species phase change (`Lee`)**: Implements volumetric evaporation and boiling coupled with multi-component species transport across solid-fluid conjugate boundaries (in contrast to free-surface VOF solvers such as `interCondensatingEvaporatingFoam` which resolve sharp liquid-gas interfaces).

---

## Benchmark Showcase: Gallium Melting (Brent et al. 1988)

The solver suite is verified against the 2D pure gallium melting benchmark (Ra ≈ 6 × 10⁵, Pr = 0.0216, Ste = 0.039). Heated from the left vertical wall (38.0°C) into an initially subcooled solid cavity (28.3°C, Tm = 29.78°C), buoyancy-driven natural convection develops a clockwise circulation roll, producing the characteristic top-dominated melt morphology.

![Gallium 2D Melting Animation](docs/images/gallium_melting_brent_comparison.gif)

![Gallium Melt Front & Liquid Fraction Comparison](docs/images/gallium_melt_front_comparison.png)

* **Melt front tracking**: The computed interface tracks digitized benchmark contours from Brent et al. (1988, Fig. 9) with mean spatial deviations of 1.7 mm at 2 min, 1.9 mm at 6 min, 2.4 mm at 10 min, and 3.5 mm at 17 min (3.89% of cavity width). These deviations are comparable to the front's diffuse numerical thickness ($\Delta T_m = 0.5$ K, cell width $\Delta x \approx 2.1$ mm).
* **Liquid fraction progression**: Simulated $f_L$ follows the Ye and Arıcı (2023) benchmark correlation (0.5736 vs 0.5703 at 17 min), while sitting above Brent et al.'s 1988 coarse numerical prediction (0.5334).
* **Energy conservation**: `EHCDarcy` achieves discrete energy closure of -0.13% over the full 17-minute transient ($42 \times 32$ grid, $n_{\text{outer}} = 8$, $\Delta t = 0.1$ s) evaluated via right-endpoint Riemann flux summation. For fast 2-minute screening runs, $n_{\text{outer}} = 3$ achieves -0.38% closure in 15.6 s CPU time.

---

## Limitations

* **Validation scope**: Benchmark verification is currently limited to 2D laminar natural convection (Rung 1) and 2D laminar melting (Rung 2). High Rayleigh number turbulent regimes are not modeled.
* **Fluid-region PCM property constraint**: In fluid regions with solid-liquid PCM, thermophysical properties must satisfy $\rho_s = \rho_l = \rho_{\text{thermo}}$ and $C_{p,s} = C_{p,l} = C_{p,\text{thermo}}$ (constant-property Boussinesq approximation).
* **Lee model kinetics**: Models volumetric boiling and condensation driven by deviation from total-pressure saturation temperature $T_{\text{sat}}(p)$; it does not model sub-boiling partial-pressure surface evaporation or free-surface liquid-gas interfaces.
* **Parallel execution**: Multi-region parallel decomposition (`decomposePar`) is currently untested; all verification runs have been performed in serial execution.
* **Inter-region mass transfer**: Phase change mass transfer is confined within each region; mass does not cross conjugate solid-fluid boundaries.
* **Diffuse interface**: Phase boundaries are represented across a numerical transition window ($\Delta T_m$), not as discontinuous sharp interfaces.

---

## Requirements and Environment

* **OpenFOAM**: OpenFOAM-v2412 (tested on Ubuntu 24.04 LTS x86_64, GCC 13.3).
* **Python**: Python 3.8+ with NumPy, Matplotlib, and Pillow (PIL) (for automated verification suites and report generation).

---

## Build Instructions

Source your OpenFOAM installation and compile with `./Allwmake`:

```bash
# Source OpenFOAM environment (path depends on your installation)
source /usr/lib/openfoam/openfoam2412/etc/bashrc

# Compile solvers and utilities
./Allwmake
```

Compiled binaries (`phaseChangeMultiRegionFoam`, `phaseChangeMultiRegionSimpleFoam`, and `writePCDict`) are installed into `$FOAM_USER_APPBIN`. To clean build artifacts, run `./Allwclean`.

### Running with Docker

A container image pre-equipped with OpenFOAM-v2412 and compiled `phaseChangeFoam` binaries is available via GitHub Container Registry:

```bash
# Pull and start an interactive simulation shell mounting current directory
docker run --rm -it -v "$PWD":/home/openfoam/run -w /home/openfoam/run \
    ghcr.io/chinmay-s-patil/phasechangefoam:latest
```

---

## Quick Start and Verification

### Running a Tutorial Case

Run the organic PCM melting tutorial (`pcm_ehc_savE_OM37`) with subcooling hysteresis:

```bash
cd tutorials/pcm_ehc_savE_OM37
blockMesh -region pcm
phaseChangeMultiRegionFoam
```

* **Expected runtime**: ~8 seconds on a single CPU core.
* **Post-processing**: Open the case in ParaView (`paraFoam -region pcm` or touch `case.foam`) and view fields `pcm/phaseFraction` and `pcm/T` to observe the advancing melt front and phase state transitions.

### Running Test Suites

```bash
# 1. Master verification suite (13 suites, 60+ quantitative checks, ~90 s)
python3 verification/runAll.py

# 2. Validation ladder fast gate (Rung 1 natural convection, ~2 s)
python3 validation/runAll.py --fast

# 3. Full validation ladder (Rungs 1 and 2, ~37 s)
python3 validation/runAll.py
```

---

## Documentation Index

Detailed reference documentation is versioned in the `docs/` directory:

| Document | Description |
| :--- | :--- |
| [`docs/models.md`](docs/models.md) | Governing equations, thermodynamic consistency, species evaporation, and PCM formulations |
| [`docs/dictionary.md`](docs/dictionary.md) | Keyword reference tables, minimal templates per model, and solver controls |
| [`docs/tutorials.md`](docs/tutorials.md) | Setup walkthroughs for boiling and PCM melting cases with manufacturer datasheets |
| [`docs/validation.md`](docs/validation.md) | Verification ladder, Brent gallium melting benchmark, and Roache GCI spatial study |
| [`docs/development.md`](docs/development.md) | Code architecture, PIMPLE execution sequence, Picard loop, and adding new models |

---

## References

To cite the `phaseChangeFoam` software, see [`CITATION.cff`](CITATION.cff).

### Physical and Numerical Model References

* Bonacina, C., Comini, G., Fasano, A., and Primicerio, M. (1973). "Numerical solution of phase-change problems", *Int. J. Heat Mass Transfer*, 16(10), pp. 1825–1832. [DOI: 10.1016/0017-9310(73)90202-0](https://doi.org/10.1016/0017-9310(73)90202-0).
* Brent, A. D., Voller, V. R., and Reid, K. J. (1988). "Enthalpy-Porosity Technique for Modeling Convection-Diffusion Phase Change: Application to the Melting of a Pure Metal", *Numerical Heat Transfer*, 13(3), pp. 297–318. [DOI: 10.1080/10407788808913615](https://doi.org/10.1080/10407788808913615).
* Carman, P. C. (1937). "Fluid flow through granular beds", *Trans. Inst. Chem. Eng.*, 15, pp. 150–166.
* Hashemi, H. T., and Sliepcevich, C. M. (1967). "A Numerical Method for Solving Two-Dimensional Problems of Heat Conduction with Change of Phase", *Chem. Eng. Prog. Symp. Ser.*, 63(79), pp. 34–41.
* Lee, W. H. (1980). "A Pressure Iteration Scheme for Two-Phase Flow Modeling", *Multiphase Transport: Fundamentals, Reactor Safety, Applications*, 1, pp. 407–432.
* Voller, V. R., and Prakash, C. (1987). "A fixed grid numerical modelling methodology for convection-diffusion mushy region phase-change problems", *Int. J. Heat Mass Transfer*, 30(8), pp. 1709–1719. [DOI: 10.1016/0017-9310(87)90317-6](https://doi.org/10.1016/0017-9310(87)90317-6).

### Benchmark and Verification References

* de Vahl Davis, G. (1983). "Natural convection of air in a square cavity: a bench mark numerical solution", *Int. J. Numer. Meth. Fluids*, 3(3), pp. 249–264. [DOI: 10.1002/fld.1650030305](https://doi.org/10.1002/fld.1650030305).
* Gau, C., and Viskanta, R. (1986). "Melting and solidification of a pure metal on a vertical wall", *J. Heat Transfer*, 108(1), pp. 174–181. [DOI: 10.1115/1.3246884](https://doi.org/10.1115/1.3246884).
* Hannoun, N., Alexiades, V., and Mai, T. Z. (2003). "Resolving the Controversy over Tin and Gallium Melting in a Rectangular Cavity Heated from the Side", *Numerical Heat Transfer, Part B: Fundamentals*, 44(3), pp. 253–276. [DOI: 10.1080/10407790390226974](https://doi.org/10.1080/10407790390226974).
* Roache, P. J. (1998). *Verification and Validation in Computational Science and Engineering*, Hermosa Publishers, Albuquerque, NM.
* Ye, H., and Arıcı, M. (2023). "Redefined-interface-error-based 2D verification and validation study on natural convection melting: Benchmark with pure Gallium", *Int. Commun. Heat Mass Transfer*, 147, 106952. [DOI: 10.1016/j.icheatmasstransfer.2023.106952](https://doi.org/10.1016/j.icheatmasstransfer.2023.106952).

---

## License and Trademark Notice

`phaseChangeFoam` is free software distributed under the GNU General Public License v3 ([`LICENSE`](LICENSE) or [GPL-3.0](https://www.gnu.org/licenses/gpl-3.0.html)).

*OPENFOAM® is a registered trade mark of OpenCFD Limited, producer and distributor of the OpenFOAM software via www.openfoam.com. This offering is not approved or endorsed by OpenCFD Limited.*
