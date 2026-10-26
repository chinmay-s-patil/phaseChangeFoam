# PCM Tutorial 1: Effective Heat Capacity (EHC) Model with savE® OM37

This tutorial demonstrates melting and thermal energy storage in **savE® OM37** (organic phase change material) using the **Effective Heat Capacity (EHC)** model.

## Material Data (savE® OM37 Datasheet)

- **Melting Transition**: $32^\circ\text{C}$ to $42^\circ\text{C}$ ($305.15\text{ K}$ to $315.15\text{ K}$)
- **Freezing Transition**: $42^\circ\text{C}$ to $32^\circ\text{C}$ ($315.15\text{ K}$ to $305.15\text{ K}$)
- **Latent Heat**: $186\text{ kJ/kg}$ (melting), $196\text{ kJ/kg}$ (freezing)
- **Density**: $\rho_{solid} = 973\text{ kg/m}^3$ (@ 30 °C), $\rho_{liquid} = 860\text{ kg/m}^3$ (@ 47 °C)
- **Specific Heat**: $C_{p,solid} = 2550\text{ J/(kg K)}$, $C_{p,liquid} = 2630\text{ J/(kg K)}$
- **Thermal Conductivity**: $k_{solid} = 0.16\text{ W/(m K)}$, $k_{liquid} = 0.13\text{ W/(m K)}$

## Solver Configuration

- **Solver**: `phaseChangeMultiRegionFoam` or `phaseChangeMultiRegionSimpleFoam`
- **Model**: `phaseChangeMode EHC;`
- **Hysteresis**: `active true;`

## How to Run

1. Generate your mesh in `constant/pcm/polyMesh` (e.g. using `blockMesh`).
2. Run the solver:
   ```bash
   phaseChangeMultiRegionFoam
   ```
