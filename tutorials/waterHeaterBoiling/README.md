# Water Heater Boiling Case Tutorial

This tutorial demonstrates a multi-region conjugate heat transfer (CHT) boiling simulation using `phaseChangeMultiRegionFoam` (or steady-state `phaseChangeMultiRegionSimpleFoam`).

## Case Overview

- **Fluid Region (`fluid`)**: Liquid water initially at $T = 365\text{ K}$ ($91.85^\circ\text{C}$). As heat transfers into the fluid from the heater, water reaches the saturation temperature ($T_{sat} = 373.15\text{ K}$) and evaporates into steam (`vapor`) using the **Lee phase change model**.
- **Solid Region (`heater`)**: Copper heater element maintained at $T = 400\text{ K}$ ($126.85^\circ\text{C}$) at the bottom boundary.

## Key Files & Configuration

- **`constant/fluid/phaseChangeDict`**:
  ```cpp
  active          true;
  phaseChange
  {
      active          true;
      type            Lee;          // Lee evaporation/condensation model
      liquid          liquid;       // Liquid species name
      vapor           vapor;        // Vapor species name
      C_evap          10.0;         // Evaporation relaxation rate [1/s]
      C_cond          10.0;         // Condensation relaxation rate [1/s]
      latentHeat      2.26e6;       // Latent heat of vaporization [J/kg]
      Tsat            373.15;       // Saturation temperature [K]
      pRef            101325;       // Reference pressure [Pa]
  }
  ```

- **`constant/fluid/thermophysicalProperties`**: Multi-component mixture defined with species `liquid` ($C_p = 4182\text{ J/(kg K)}$) and `vapor` ($C_p = 2080\text{ J/(kg K)}$, $H_f = 2.26\times 10^6\text{ J/kg}$).
- **`constant/heater/thermophysicalProperties`**: Copper solid properties ($k = 385\text{ W/(m K)}$, $\rho = 8960\text{ kg/m}^3$, $C_p = 385\text{ J/(kg K)}$).

## How to Run

1. Generate or copy your mesh into `constant/fluid/polyMesh` and `constant/heater/polyMesh` (e.g. using `blockMesh` / `splitMeshByTopology`).
2. Run the solver:
   ```bash
   phaseChangeMultiRegionFoam
   ```
   or for steady-state estimation:
   ```bash
   phaseChangeMultiRegionSimpleFoam
   ```
