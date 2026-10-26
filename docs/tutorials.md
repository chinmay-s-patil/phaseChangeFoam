# Tutorials & Case Setup

Case walkthroughs, boundary conditions, material datasheets, and execution procedures for `phaseChangeFoam`.

> **Note**: Tutorial cases are intended as demonstration workflows illustrating case configuration, meshing, and model options. They are not formally validated against experimental measurements. For quantitative verification against literature benchmarks, see [`docs/validation.md`](validation.md).

---

## Overview of Tutorial Cases

| Tutorial Case | Directory | Primary Physics | Solvers | Key Models Used |
| :--- | :--- | :--- | :--- | :--- |
| Water Heater Boiling | [`tutorials/waterHeaterBoiling/`](../tutorials/waterHeaterBoiling/) | Multi-region CHT, boiling, evaporation | `phaseChangeMultiRegionFoam` | `Lee`, Multi-component species |
| savE OM37 Melting | [`tutorials/pcm_ehc_savE_OM37/`](../tutorials/pcm_ehc_savE_OM37/) | Organic PCM thermal storage, hysteresis | `phaseChangeMultiRegionFoam` | `EHC`, Thermal hysteresis |
| savE HS36 Mushy Flow | [`tutorials/pcm_ep_savE_HS36/`](../tutorials/pcm_ep_savE_HS36/) | Inorganic PCM melting, Darcy drag | `phaseChangeMultiRegionFoam` | `enthalpyPorosity`, Darcy drag |

---

## Tutorial 1: Water Heater Boiling (`waterHeaterBoiling`)

* **Case directory**: [`tutorials/waterHeaterBoiling/`](../tutorials/waterHeaterBoiling/)
* **Estimated runtime**: ~25 seconds on a single core.

### Physical Setup
Simulates conjugate heat transfer and boiling in a two-region domain:
* **Solid region (`heater`)**: Copper block with bottom boundary maintained at $T = 400$ K (126.85°C).
* **Fluid region (`fluid`)**: Water initially subcooled at $T = 365$ K (91.85°C). Across the coupled interface, heat conducts into the fluid. Near-wall water reaches $T_{\text{sat}} = 373.15$ K and evaporates into vapor via the volumetric Lee model.

### Thermophysical Setup
Configured in `constant/fluid/thermophysicalProperties` as a multi-component reacting mixture with two species:
* `liquid`: Liquid water ($C_p = 4182$ J/(kg K), molecular weight $W = 18.015$ kg/kmol).
* `vapor`: Water vapor ($C_p = 2080$ J/(kg K), molecular weight $W = 18.015$ kg/kmol, $H_f = 2.26 \times 10^6$ J/kg).

### Running the Case
```bash
cd tutorials/waterHeaterBoiling
blockMesh
splitMeshRegions -cellZones -overwrite
phaseChangeMultiRegionFoam
```

### Post-Processing & Expected Results
Open in ParaView (`paraFoam -touch` and open the `.foam` file):
* Inspect `fluid/vapor` (mass fraction): vapor accumulates along the heated interface and diffuses into the headspace.
* Inspect `fluid/T` and `heater/T`: verify continuity of temperature and heat flux across the conjugate boundary.

---

## Tutorial 2: Organic PCM Melting with Hysteresis (`pcm_ehc_savE_OM37`)

* **Case directory**: [`tutorials/pcm_ehc_savE_OM37/`](../tutorials/pcm_ehc_savE_OM37/)
* **Estimated runtime**: ~8 seconds on a single core.

### Physical Setup & Material Datasheet
Simulates thermal storage in a bio-based organic phase change material (SavE® OM37, Pluss Advanced Technologies).
* Initial solid temperature: $T_{\text{init}} = 300$ K (26.85°C).
* Left boundary: Heated at $T_h = 320$ K (46.85°C).
* Top, bottom, right boundaries: Adiabatic.

| Property (SavE® OM37) | Value | Unit | Datasheet Source |
| :--- | :--- | :--- | :--- |
| Melting nominal temperature | 37 (310.15 K) | °C | Pluss Advanced Technologies |
| Freezing nominal temperature | 36 (309.15 K) | °C | Pluss Advanced Technologies |
| Melting window ($T_{\text{sol,m}} \to T_{\text{liq,m}}$) | 309.65 – 310.65 (36.5 – 37.5°C) | K | Nominal ±0.5 K window |
| Freezing window ($T_{\text{sol,f}} \to T_{\text{liq,f}}$) | 308.65 – 309.65 (35.5 – 36.5°C) | K | Subcooled by 1.0 K |
| Latent heat of fusion ($L$) | 210000 | J/kg | Pluss Advanced Technologies |
| Liquid density | 820 | kg/m3 | Pluss Advanced Technologies |
| Solid density | 910 | kg/m3 | Pluss Advanced Technologies |
| Specific heat capacity ($C_p$) | 2200 | J/(kg K) | Pluss Advanced Technologies |
| Thermal conductivity ($k$) | 0.15 | W/(m K) | Pluss Advanced Technologies |

### Physical Meaning of Hysteresis in OM37
SavE® OM37 exhibits subcooling during thermal cycling. During heating, the material remains solid until $36.5^\circ\text{C}$ and finishes melting at $37.5^\circ\text{C}$. During cooling, the liquid does not freeze at $37^\circ\text{C}$; instead, it supercools to $36.5^\circ\text{C}$ before latent heat release begins, completing solidification at $35.5^\circ\text{C}$. The model tracks this path dependence using the cell-wise `phaseState` variable.

### Running the Case
```bash
cd tutorials/pcm_ehc_savE_OM37
blockMesh -region pcm
phaseChangeMultiRegionFoam
```

### Post-Processing & Expected Results
Open the case in ParaView (`paraFoam -region pcm`):
* Inspect `phaseFraction`: a planar melt front advances from the left heated wall toward the right.
* Inspect `phaseState`: confirms the transition sequence from solid (`0`) through melting (`1`) to liquid (`2`).
* Inspect `CpEff`: shows a localized peak in apparent heat capacity inside the mushy zone ($C_{p,\text{eff}} \approx 210000 / 1.0 + 2200 \approx 212200$ J/(kg K)).

---

## Tutorial 3: Salt Hydrate Mushy Zone Flow (`pcm_ep_savE_HS36`)

* **Case directory**: [`tutorials/pcm_ep_savE_HS36/`](../tutorials/pcm_ep_savE_HS36/)
* **Estimated runtime**: ~12 seconds on a single core.

### Physical Setup & Material Datasheet
Simulates an inorganic salt hydrate PCM (SavE® HS36, Pluss Advanced Technologies) modeled in a fluid region using `enthalpyPorosity` with Darcy mushy-zone drag.

| Property (SavE® HS36) | Value | Unit | Datasheet Source |
| :--- | :--- | :--- | :--- |
| Nominal melting temperature | 36 (309.15 K) | °C | Pluss Advanced Technologies |
| Melting window ($T_{\text{solidus}} \to T_{\text{liquidus}}$) | 308.65 – 309.65 (35.5 – 36.5°C) | K | Nominal ±0.5 K window |
| Latent heat of fusion ($L$) | 270000 | J/kg | Pluss Advanced Technologies |
| Reference density | 1500 | kg/m3 | Pluss Advanced Technologies |
| Specific heat capacity ($C_p$) | 2300 | J/(kg K) | Pluss Advanced Technologies |
| Thermal conductivity ($k$) | 0.54 | W/(m K) | Pluss Advanced Technologies |
| Darcy coefficient ($C_u$) | 1.0e5 | kg/(m3 s) | Carman-Kozeny parameter |

### Running the Case
```bash
cd tutorials/pcm_ep_savE_HS36
blockMesh -region pcm
phaseChangeMultiRegionFoam
```

### Post-Processing & Expected Results
Open in ParaView (`paraFoam -region pcm`):
* Inspect `phaseFraction`: observe the advancing melting interface.
* Inspect `U`: velocity is damped to zero in the solid region by the Darcy drag sink $\mathbf{S}_u$.
