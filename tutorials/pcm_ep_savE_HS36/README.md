# PCM Tutorial 2: Enthalpy-Porosity (EP) Model with savE® HS36

This tutorial demonstrates melting and phase change in **savE® HS36** (inorganic salt hydrate PCM) using the **Enthalpy-Porosity (EP)** model with mushy-zone flow resistance formulation.

## Material Data (savE® HS36 Datasheet)

- **Melting Transition**: $30^\circ\text{C}$ to $40^\circ\text{C}$ ($303.15\text{ K}$ to $313.15\text{ K}$)
- **Freezing Transition**: $40^\circ\text{C}$ to $30^\circ\text{C}$ ($313.15\text{ K}$ to $303.15\text{ K}$)
- **Latent Heat**: $163\text{ kJ/kg}$ (melting and freezing)
- **Density**: $\rho_{solid} = 1967\text{ kg/m}^3$ (@ 20 °C), $\rho_{liquid} = 1850\text{ kg/m}^3$ (@ 45 °C)
- **Specific Heat**: $C_{p,solid} = 1980\text{ J/(kg K)}$, $C_{p,liquid} = 2320\text{ J/(kg K)}$
- **Thermal Conductivity**: $k_{solid} = 0.50\text{ W/(m K)}$, $k_{liquid} = 0.47\text{ W/(m K)}$

## Solver Configuration

- **Solver**: `phaseChangeMultiRegionFoam` or `phaseChangeMultiRegionSimpleFoam`
- **Model**: `phaseChangeMode enthalpyPorosity;`
- **Mushy Zone Parameters**: `A_cu 1.0e5; eps 1.0e-3;`

## How to Run

1. Generate your mesh in `constant/pcm/polyMesh` (e.g. using `blockMesh`).
2. Run the solver:
   ```bash
   phaseChangeMultiRegionFoam
   ```
