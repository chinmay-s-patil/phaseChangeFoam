# Validation Suite Overview

This directory contains a progressive five-rung assessment suite covering natural convection, cyclic thermal hysteresis state-logic verification, liquid-metal melting benchmarks, and preliminary comparison with physical experimental data.

For the complete verification and validation methodology, quantitative results, and experimental validation status, see [`docs/validation.md`](../docs/validation.md).

---

## Validation Hierarchy

```
[ Rung 5: 3D Pure Gallium Cavity Melting ]
Gau & Viskanta (1986)
Physical experimental validation screening with no-slip sidewall braking
Transient solver: phaseChangeMultiRegionFoam
               ▲
               │  (Physical experimental validation)
               │
[ Rung 4: 2D Gallium Melting Benchmark ]
Brent, Voller, Reid (1988) / Ye & Arıcı (2023)
Coupled phase change, Darcy mushy drag, natural convection (Ra ≈ 6.06e5, Pr = 0.0216, Ste = 0.039)
Transient solver: phaseChangeMultiRegionFoam
               ▲
               │  (Coupled convection-diffusion phase change)
               │
[ Rung 3: 2D Pure Tin Melting ]
Bertrand et al. (1999) / Hannoun et al. (2003)
Liquid metal melting benchmark (Ra = 2.50e5, Pr = 0.02, Ste = 0.01)
Transient solver: phaseChangeMultiRegionFoam
               ▲
               │  (Low-Prandtl number liquid metal regime)
               │
[ Rung 2: EHC Cyclic Thermal Hysteresis ]
Direction-dependent melting and freezing loop with 5 K hysteresis window
State-variable subcooling tracking and discrete energy conservation
Transient solver: phaseChangeMultiRegionFoam
               ▲
               │  (Hysteresis path state logic)
               │
[ Rung 1: Natural Convection Cavity ]
de Vahl Davis (1983)
Pure buoyancy circulation (Ra = 10^5, Pr = 0.71) in square cavity
Steady-state solver: phaseChangeMultiRegionSimpleFoam
```

---

## Benchmarks Overview

### Rung 1: de Vahl Davis (1983) Natural Convection Cavity
- **Directory**: [`deVahlDavis2D/`](deVahlDavis2D/)
- **Solver**: `phaseChangeMultiRegionSimpleFoam` (steady-state Boussinesq)
- **Physics**: Laminar natural convection in a square cavity ($L = 0.01\text{ m}$) at $\text{Ra} = 10^5, \text{Pr} = 0.71$.
- **Headline Result**: Hot-wall mean Nusselt number matches benchmark within 0.70% on $50 \times 50$ ($\overline{\text{Nu}} = 4.551$ vs $4.519$) and 0.64% on $80 \times 80$.
- **Runtime**: $\approx 2\text{ s}$ (500 iterations).

### Rung 2: EHC Cyclic Thermal Hysteresis
- **Directory**: [`pcmHysteresis/`](pcmHysteresis/)
- **Solver**: `phaseChangeMultiRegionFoam`
- **Physics**: Direction-dependent cyclic melting ($[308, 312]\text{ K}$) and freezing ($[303, 307]\text{ K}$) loop with state-variable subcooling retention.
- **Headline Result**: Complete phase cycle recovery with subcooling delay (100% molten $\to$ subcooling retention $\to$ 100% solid crystallization).
- **Runtime**: $\approx 2\text{ s}$.

### Rung 3: Bertrand / Hannoun Pure Tin Melting
- **Directory**: [`tin2D/`](tin2D/)
- **Solver**: `phaseChangeMultiRegionFoam`
- **Physics**: Low-Prandtl liquid metal melting in a square cavity ($\text{Ra} = 2.50 \times 10^5, \text{Pr} = 0.02, \text{Ste} = 0.01$).
- **Headline Result**: Early conduction front ($t = 100\text{ s}$) matches the 1D Stefan limit within $-1.34\%$ ($f_L = 0.0881$ vs $0.0893$). Full late-time convective benchmark ($t = 2500\text{ s}$) is in progress (literature targets $f_L \approx 0.54\text{--}0.60$ depend on resolved $\text{Pr}$).
- **Runtime**: $\approx 15\text{ s}$.

### Rung 4: Brent (1988) / Ye & Arıcı (2023) 2D Gallium Melting
- **Directory**: [`gallium2D/`](gallium2D/)
- **Solver**: `phaseChangeMultiRegionFoam` (transient PIMPLE)
- **Physics**: Coupled convection and melting in a rectangular gallium cavity ($W = 8.89\text{ cm}, H = 6.35\text{ cm}$, $\text{Ra} \approx 6.06 \times 10^5$, $\text{Pr} = 0.0216$).
- **Headline Result**: Late-time liquid fraction ($f_L = 0.5736$ at $17\text{ min}$) agrees within $+0.58\%$ of the fine-grid numerical correlation of Ye & Arıcı ($0.5703$); 2-minute gate achieves discrete energy closure within $-0.34\%$.
- **Runtime**: $\approx 35\text{ s}$ for 2-min gate; $\approx 4\text{--}5\text{ min}$ for full 17-min transient.

### Rung 5: Gau & Viskanta (1986) 3D Gallium Cavity Validation Screening
- **Directory**: [`gallium3D_experimental/`](gallium3D_experimental/)
- **Solver**: `phaseChangeMultiRegionFoam` (3D)
- **Physics**: 3D rectangular cavity ($8.89\text{ cm} \times 6.35\text{ cm} \times 3.81\text{ cm}$) with physical no-slip inspection windows exerting sidewall viscous drag.
- **Headline Result**: Preliminary screening run demonstrates active 3D boundary-layer braking; formal ASME V&V 20 experimental validation is underway.
- **Runtime**: $\approx 2.5\text{ min}$.

---

## Reference Data & Literature Sources
- **Brent, Voller, Reid (1988)**: "Enthalpy-Porosity Technique for Modeling Convection-Diffusion Phase Change: Application to the Melting of a Pure Metal", *Numerical Heat Transfer*, 13(3), pp. 297–318. DOI: [10.1080/10407788808913615](https://doi.org/10.1080/10407788808913615).
- **de Vahl Davis (1983)**: "Natural convection of air in a square cavity: a bench mark numerical solution", *Int. J. Numer. Meth. Fluids*, 3, pp. 249–264. DOI: [10.1002/fld.1650030305](https://doi.org/10.1002/fld.1650030305).
- **Hannoun, Alexiades, Mai (2003)**: "Resolving the Controversy over Tin and Gallium Melting in a Rectangular Cavity Heated from the Side", *Numerical Heat Transfer, Part B: Fundamentals*, 44(3), pp. 253–276. DOI: [10.1080/10407790390226974](https://doi.org/10.1080/10407790390226974).
- **Ye & Arıcı (2023a)**: "3D validation, 2D feasibility, corrected and developed correlations for pure solid-gallium phase change modeling by enthalpy-porosity methodology", *International Communications in Heat and Mass Transfer*, 144, 106780. DOI: [10.1016/j.icheatmasstransfer.2023.106780](https://doi.org/10.1016/j.icheatmasstransfer.2023.106780). *(Source of Eq. 6 liquid fraction benchmark correlation)*.
- **Ye & Arıcı (2023b)**: "Redefined interface error, 2D verification and validation for pure solid-gallium phase change modeling by enthalpy-porosity methodology", *International Communications in Heat and Mass Transfer*, 147, 106952. DOI: [10.1016/j.icheatmasstransfer.2023.106952](https://doi.org/10.1016/j.icheatmasstransfer.2023.106952). *(Source of 2D grid-converged interface benchmarks)*.
- **Gau & Viskanta (1986)**: "Melting and Solidification of a Metal System—I. Melting from a Vertical Wall", *Int. J. Heat Mass Transfer*, 29(12), pp. 1951–1963. DOI: [10.1016/0017-9310(86)90011-2](https://doi.org/10.1016/0017-9310(86)90011-2). Traced midplane targets from Fig. 5 ($t = 30, 60, 120, 360, 600, 1020\text{ s}$) carry an estimated digitization uncertainty of $\pm 1.0\text{ mm}$.
- **`gallium2D/reference_data/`**: Digitized melt front coordinates ($t = 2, 6, 10, 17\text{ min}$) traced from Brent et al. (1988) Fig. 9 with $\pm 1\text{ mm}$ horizontal front position uncertainty.

---

## Running the Validation Ladder

Run the validation ladder from the repository root:

```bash
# Run the complete ladder (All 5 rungs, ~3.5 min):
python3 validation/runAll.py

# Run fast gate (Rungs 1 & 2: de Vahl Davis & EHC hysteresis, ~5 s):
python3 validation/runAll.py --fast

# Run individual benchmarks:
python3 validation/deVahlDavis2D/runTest.py
python3 validation/pcmHysteresis/runTest.py
python3 validation/tin2D/runTest.py
python3 validation/gallium2D/runTest.py --endTime 120 --plot
python3 validation/gallium3D_experimental/runTest.py

# Clean all generated cases and simulation runs:
./validation/Allclean
```

