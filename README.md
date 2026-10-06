# `phaseChangeFoam`: Multi-Region Phase Change & Conjugate Heat Transfer Suite

`phaseChangeFoam` is a comprehensive OpenFOAM solver suite and model library designed for multi-region conjugate heat transfer (CHT), species phase change (evaporative/boiling volumetric phase change), and solid-liquid thermal phase change (solidification, melting, thermal energy storage, and mushy zone flows).

---

## Table of Contents

1. [Overview & Solvers](#overview--solvers)
2. [Model Architecture](#model-architecture)
3. [Fluid Region Phase Change Models](#fluid-region-phase-change-models)
   - [Lee Model (`Lee`)](#1-lee-model-lee)
   - [Constant Source Model (`constantSource`)](#2-constant-source-model-constantsource)
   - [Inactive Model (`none`)](#3-inactive-model-none)
4. [Solid Region & PCM Phase Change Models](#solid-region--pcm-phase-change-models)
   - [Effective Heat Capacity Model (`EHC`)](#1-effective-heat-capacity-model-ehc)
   - [Enthalpy-Porosity Model (`enthalpyPorosity`)](#2-enthalpy-porosity-model-enthalpyporosity)
   - [Inactive Solid Model (`none`)](#3-inactive-solid-model-none)
5. [Complete Dictionary Configuration Reference (`phaseChangeDict`)](#complete-dictionary-configuration-reference-phasechangedict)
6. [C++ Code Snippets & API Usage Guide](#c-code-snippets--api-usage-guide)
   - [Fluid Region Integration](#fluid-region-c-integration)
   - [Solid Region Integration](#solid-region-c-integration)
7. [Physics Caveats & Numerical Considerations](#physics-caveats--numerical-considerations)
8. [Utilities, Tutorials, & Verification Suite](#utilities-tutorials--verification-suite)

---

## Overview & Solvers

The suite provides primary multi-region solvers and companion utility binaries:

* **[`phaseChangeMultiRegionFoam`](file:///home/lavender/OpenFoamUbu/solvers/phaseChangeFoam/src/solvers/phaseChangeMultiRegionFoam/phaseChangeMultiRegionFoam.C)**: Transient solver for multi-region conjugate heat transfer with species fluid phase change (Lee model) and solid region phase change (EHC / Enthalpy-Porosity).
* **[`phaseChangeMultiRegionSimpleFoam`](file:///home/lavender/OpenFoamUbu/solvers/phaseChangeFoam/src/solvers/phaseChangeMultiRegionSimpleFoam)**: Steady-state solver for multi-region conjugate heat transfer and phase change.
* **[`writePCDict`](file:///home/lavender/OpenFoamUbu/solvers/phaseChangeFoam/src/utilities/writePCDict/writePCDict.C)**: Utility to automatically generate fully documented template `phaseChangeDict` files for case setup.

---

## Model Architecture

The framework relies on runtime selection tables for both fluid and solid region phase change models:

```
                                      +------------------------------------+
                                      |         phaseChangeDict            |
                                      +------------------------------------+
                                                        |
                         +------------------------------+------------------------------+
                         |                                                             |
                         v                                                             v
        +----------------------------------+                         +----------------------------------+
        |   fluidPhaseChangeModel (Base)   |                         |     phaseChangeModel (Base)      |
        +----------------------------------+                         +----------------------------------+
                         |                                                             |
    +--------------------+--------------------+                   +--------------------+--------------------+
    |                    |                    |                   |                    |                    |
    v                    v                    v                   v                    v                    v
+-------+        +----------------+       +------+            +-------+        +------------------+     +------+
|  Lee  |        | constantSource |       | none |            |  EHC  |        | enthalpyPorosity |     | none |
+-------+        +----------------+       +------+            +-------+        +------------------+     +------+
```

* **Fluid Models** derive from [`fluidPhaseChangeModel`](file:///home/lavender/OpenFoamUbu/solvers/phaseChangeFoam/src/solvers/phaseChangeMultiRegionFoam/fluid/fluidPhaseChangeModel.H) and feed mass sources to continuity/pressure (`pEqn`), species mass fraction sources to `YEqn`, and latent heat sources to energy equations (`EEqn`).
* **Solid/PCM Models** derive from [`phaseChangeModel`](file:///home/lavender/OpenFoamUbu/solvers/phaseChangeFoam/src/solvers/phaseChangeMultiRegionFoam/solid/phaseChangeModel.H) and supply path-integrated heat capacities $C_{p,eff}$, exact Newton secant linearisation $S_p$, latent energy sources $S_u$, effective thermal conductivities $k_{eff}$, densities $\rho_{eff}$, and Darcy momentum drag forces $S_{u,drag}$.

---

## Fluid Region Phase Change Models

Activated in fluid regions via `constant/<region>/phaseChangeDict`.

### 1. Lee Model (`Lee`)

Class: [`leeFluidPhaseChangeModel`](file:///home/lavender/OpenFoamUbu/solvers/phaseChangeFoam/src/solvers/phaseChangeMultiRegionFoam/fluid/leeFluidPhaseChangeModel.H)

The volumetric Lee model calculates species mass transfer rate $\dot{m}$ ($\text{kg}/(\text{m}^3\,\text{s})$) between liquid and vapor species based on local temperature relative to saturation temperature $T_{sat}$.

#### Mathematical Formulation & Explicit Capping
To guarantee PIMPLE outer-corrector invariance across $nOuterCorrectors$, evaluation uses old-time state $(T_0, p_0, \rho_0, Y_{l,0}, Y_{v,0})$ during outer corrector step `oCorr == 0`:

* **Pressure-Dependent Saturation Temperature $T_{sat}(p)$** (when `enableTsatP true`):
  $$\frac{1}{T_{sat}(p)} = \frac{1}{T_{sat,0}} - \frac{R_v}{L} \ln\left(\frac{p_0}{p_{ref}}\right), \quad R_v = \frac{R_{universal}}{W_v}$$

* **Evaporation Rate ($T_0 > T_{sat}$)**:
  $$\dot{m}_{evap} = \min\left( C_{evap} \rho_0 Y_{l,0} \frac{T_0 - T_{sat}}{T_{sat}},\, \frac{\rho_0 Y_{l,0}}{\Delta t},\, \frac{\rho_0 C_p (T_0 - T_{sat})}{L \Delta t} \right)$$

* **Condensation Rate ($T_0 < T_{sat}$)**:
  $$\dot{m}_{cond} = - \min\left( C_{cond} \rho_0 Y_{v,0} \frac{T_{sat} - T_0}{T_{sat}},\, \frac{\rho_0 Y_{v,0}}{\Delta t},\, \frac{\rho_0 C_p (T_{sat} - T_0)}{L \Delta t} \right)$$

* **Energy Source Term**:
  $$\text{energySource} = - \dot{m} \cdot L \quad [\text{W/m}^3]$$

#### Lee Model Configuration Options
```foam
phaseChange
{
    active          true;
    type            Lee;            // Selection keyword
    liquid          H2O_l;          // Name of liquid specie in thermo composition
    vapor           H2O_v;          // Name of vapor specie in thermo composition
    C_evap          0.1;            // Evaporation relaxation frequency coefficient [1/s]
    C_cond          0.1;            // Condensation relaxation frequency coefficient [1/s]
    latentHeat      2.26e6;         // Latent heat of vaporization L [J/kg]
    Tsat            373.15;         // Reference saturation temperature Tsat [K]
    enableTsatP     false;          // Toggle pressure-dependent Tsat(p) via Clausius-Clapeyron
    pRef            101325;         // Reference saturation pressure pRef [Pa]
}
```

---

### 2. Constant Source Model (`constantSource`)

Class: [`constantSourceFluidPhaseChangeModel`](file:///home/lavender/OpenFoamUbu/solvers/phaseChangeFoam/src/solvers/phaseChangeMultiRegionFoam/fluid/constantSourceFluidPhaseChangeModel.H)

Verification model providing user-specified constant volumetric sources for plumbing testing.

```foam
phaseChange
{
    active          true;
    type            constantSource;
    speciesName     H2O;            // Specie name receiving mass source
    massSource      0.0;            // Net mass source [kg/(m^3 s)]
    energySource    -1.0e6;         // Energy source [W/m^3]
    speciesSource   1.0e-3;         // Species mass source [kg/(m^3 s)]
}
```

---

### 3. Inactive Model (`none`)

Class: [`noneFluidPhaseChangeModel`](file:///home/lavender/OpenFoamUbu/solvers/phaseChangeFoam/src/solvers/phaseChangeMultiRegionFoam/fluid/noneFluidPhaseChangeModel.H)

Returns zero sources for mass, species, and energy equations.

```foam
phaseChange
{
    active          false;
    type            none;
}
```

---

## Solid Region & PCM Phase Change Models

Activated in solid/PCM regions via `constant/<region>/phaseChangeDict`.

### 1. Effective Heat Capacity Model (`EHC`)

Class: [`ehcPhaseChangeModel`](file:///home/lavender/OpenFoamUbu/solvers/phaseChangeFoam/src/solvers/phaseChangeMultiRegionFoam/solid/ehcPhaseChangeModel.H)

Solid-liquid phase change model utilizing path-integral Effective Heat Capacity ($C_{p,eff}$) with exact Newton/secant linearisation of latent heat, stateful hysteresis tracking, and path-integral thermophysical property evaluation.

#### Key Features:
* **Transformed Phase Fraction $\alpha \in [0, 1]$**: Smooth transition between solidus ($T_{solidus}$) and liquidus ($T_{liquidus}$).
* **Direction Restriction (`direction`)**:
  - `both`: Bidirectional melting and freezing with optional thermal hysteresis.
  - `forward`: Irreversible melting (phase fraction cannot decrease).
  - `reverse`: Irreversible freezing (phase fraction cannot increase).
* **Stateful Thermal Hysteresis (`hysteresis`)**: Tracks heating vs cooling history ($T_{reversal}$) and enforces hysteresis plateaus when temperature reverses within deadband `reversalTolerance`.
* **Path-Integral Thermophysics (`thermophysical`)**: Computes exact step sensible enthalpy change $\Delta h_{sens} / \Delta T$ to prevent numerical temperature overshoot across phase transition boundaries.
* **Exact Newton Secant Linearization**: Generates implicit diagonal coefficient $S_p = \frac{\bar{\rho} L \frac{d\alpha}{dT}}{\Delta t C_{p,thermo}}$ and explicit source $S_u = \frac{\bar{\rho} L (\alpha - \alpha_{old})}{\Delta t}$.

---

### 2. Enthalpy-Porosity Model (`enthalpyPorosity`)

Class: [`enthalpyPorosityPhaseChangeModel`](file:///home/lavender/OpenFoamUbu/solvers/phaseChangeFoam/src/solvers/phaseChangeMultiRegionFoam/solid/enthalpyPorosityPhaseChangeModel.H)

Extends the solid phase change model for fluid mushy zones by adding Darcy momentum damping force to penalize velocity in solidifying/melting regions.

#### Damping Force Formulation:
$$\vec{S}_{u,drag} = - C_u \frac{(1 - \alpha)^2}{\alpha^3 + q} \vec{U} \quad \left[\frac{\text{N}}{\text{m}^3}\right]$$

Where:
* $C_u$ (or `A_cu`, `Cu`): Mushy zone Darcy constant ($\text{kg}/(\text{m}^3\,\text{s})$), default `1e5`.
* $q$ (or `eps`): Division tolerance to prevent division by zero, default `0.001`.

---

### 3. Inactive Solid Model (`none`)

Class: [`nonePhaseChangeModel`](file:///home/lavender/OpenFoamUbu/solvers/phaseChangeFoam/src/solvers/phaseChangeMultiRegionFoam/solid/nonePhaseChangeModel.H)

Disables solid phase change calculations for pure conduction or single-phase solid regions.

---

## Complete Dictionary Configuration Reference (`phaseChangeDict`)

Below is a complete, fully documented template for `constant/<region>/phaseChangeDict` containing all available options:

```foam
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "constant/air";
    object      phaseChangeDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

active          true;

phaseChange
{
    active          true;

    // =========================================================================
    // 1. Solid / PCM Phase Change Model & Direction Controls
    // =========================================================================
    type            EHC;        // Selection mode: EHC, enthalpyPorosity, none
    direction       both;       // Options: both, forward (melting only), reverse (freezing only)

    // Forward / Melting Phase Transition
    forward
    {
        T_lowerBound    300.0;  // Solidus temperature T_solidus [K]
        T_upperBound    310.0;  // Liquidus temperature T_liquidus [K]
        latentHeat      100000; // Latent heat of fusion L [J/kg]
    }

    // Reverse / Freezing Phase Transition (Thermal Hysteresis)
    reverse
    {
        T_lowerBound    295.0;  // Freezing lower bound temperature [K]
        T_upperBound    305.0;  // Freezing upper bound temperature [K]
        latentHeat      100000; // Latent heat of freezing L [J/kg]
    }

    // Stateful Thermal Hysteresis Tracking
    hysteresis
    {
        active              false;  // Enable stateful hysteresis tracking
        reversalTolerance   1e-6;   // Temperature reversal deadband tolerance [K]
    }

    // Phase Density Model
    density
    {
        model                       thermo; // Options: thermo (constant), linear
        rhoRef                      1000.0; // Reference density [kg/m^3]
        rhoSolid                    1000.0; // Solid phase density [kg/m^3]
        rhoLiquid                   1000.0; // Liquid phase density [kg/m^3]
        allowNonConservativeDensity false;  // Allow unequal solid/liquid density
    }

    // Thermophysical Properties Mode
    thermophysical
    {
        mode            thermo; // Options: thermo (from thermo dict), custom
        CpSolid         1980.0; // Custom solid specific heat capacity Cp [J/(kg K)]
        CpLiquid        2320.0; // Custom liquid specific heat capacity Cp [J/(kg K)]
        kSolid          1.0;    // Custom solid thermal conductivity k [W/(m K)]
        kLiquid         1.0;    // Custom liquid thermal conductivity k [W/(m K)]
    }

    // Convection Suppression (Solid Regions)
    convection
    {
        suppress        true;   // Suppress fluid convection in solid regions
    }

    // Mushy Zone Porosity Drag (Enthalpy-Porosity Mode Only)
    porosity
    {
        Cu              1.0e5;  // Darcy penalty coefficient Cu [kg/(m^3 s)]
        q               1.0e-3; // Small denominator tolerance q
    }

    // Thermal Buoyancy Coefficient
    buoyancy
    {
        beta            0.0;    // Thermal expansion coefficient beta [1/K]
    }

    // =========================================================================
    // 2. Fluid Phase Change Model (Volumetric Evaporation / Condensation)
    // =========================================================================
    liquid          H2O_l;      // Liquid specie name in thermo composition
    vapor           H2O_v;      // Vapor specie name in thermo composition
    C_evap          0.1;        // Evaporation rate coefficient [1/s]
    C_cond          0.1;        // Condensation rate coefficient [1/s]
    latentHeat      2.26e6;     // Latent heat of vaporization L [J/kg]
    Tsat            373.15;     // Reference saturation temperature Tsat [K]
    enableTsatP     false;      // Enable Clausius-Clapeyron Tsat(p)
    pRef            101325;     // Reference saturation pressure pRef [Pa]
}
```

### Table of Options & Defaults

| Block / Keyword | Description | Type | Default |
| :--- | :--- | :--- | :--- |
| `phaseChange.active` | Master toggle for phase change model | `bool` | `true` |
| `phaseChange.type` | Model selection (`Lee`, `constantSource`, `EHC`, `enthalpyPorosity`, `none`) | `word` | `Lee` / `EHC` |
| `phaseChange.direction` | Direction restriction (`both`, `forward`, `reverse`) | `word` | `both` |
| `forward.T_lowerBound` | Solidus / lower bound temperature | `scalar` [K] | `300.0` |
| `forward.T_upperBound` | Liquidus / upper bound temperature | `scalar` [K] | `310.0` |
| `forward.latentHeat` | Latent heat of transition $L$ | `scalar` [J/kg] | `1.0e5` |
| `hysteresis.active` | Enables path-dependent thermal hysteresis | `bool` | `false` |
| `hysteresis.reversalTolerance` | Deadband tolerance for temperature direction reversal | `scalar` [K] | `1e-6` |
| `density.model` | Phase density interpolation mode (`thermo`, `linear`) | `word` | `thermo` |
| `thermophysical.mode` | Thermophysical mode (`thermo`, `custom`) | `word` | `thermo` |
| `porosity.Cu` | Darcy drag coefficient $C_u$ | `scalar` [kg/(m³ s)] | `1.0e5` |
| `porosity.q` | Division prevention tolerance $q$ | `scalar` | `0.001` |
| `C_evap` / `C_cond` | Lee evaporation and condensation rate coefficients | `scalar` [1/s] | `0.1` |
| `enableTsatP` | Clausius-Clapeyron pressure-dependent $T_{sat}(p)$ toggle | `bool` | `false` |

---

## C++ Code Snippets & API Usage Guide

### Fluid Region C++ Integration

In a fluid region solver (e.g. [`solveFluid.H`](file:///home/lavender/OpenFoamUbu/solvers/phaseChangeFoam/src/solvers/phaseChangeMultiRegionFoam/fluid/solveFluid.H)):

```cpp
#include "fluidPhaseChangeModel.H"

// 1. Instantiation via Runtime Selection Table
autoPtr<fluidPhaseChangeModel> phaseChange = fluidPhaseChangeModel::New
(
    mesh,
    thermo,
    U,
    phi
);

// 2. Correct state and verify mass conservation in solver loop
phaseChange->correct();
phaseChange->checkMassConservation();

// 3. Continuity / Pressure Equation (pEqn.H)
if (phaseChange->active())
{
    // Subtract mass source from continuity equation
    p_rghDDtEqn -= phaseChange->massSource();
}

// 4. Species Transport Equation (YEqn.H)
forAll(species, i)
{
    fvScalarMatrix YiEqn
    (
        fvm::ddt(rho, Yi)
      + mvConvection->fvmDiv(phi, Yi)
      - fvm::laplacian(turbulence->muEff(), Yi)
    );

    if (phaseChange->active())
    {
        // Add species source/sink term
        YiEqn -= phaseChange->speciesSource(i);
    }

    YiEqn.solve();
}

// 5. Fluid Energy Equation (EEqn.H)
fvScalarMatrix EEqn
(
    fvm::ddt(rho, he)
  + mvConvection->fvmDiv(phi, he)
  + /.../
);

if (phaseChange->active())
{
    // Subtract latent heat source term
    EEqn -= phaseChange->energySource();
}
EEqn.solve();
```

---

### Solid Region C++ Integration

In a solid region solver (e.g. [`phaseChangeMultiRegionFoam.C`](file:///home/lavender/OpenFoamUbu/solvers/phaseChangeFoam/src/solvers/phaseChangeMultiRegionFoam/phaseChangeMultiRegionFoam.C)):

```cpp
#include "phaseChangeModel.H"

// 1. Instantiation via Runtime Selection Table
autoPtr<phaseChangeModel> phaseChange = phaseChangeModel::New
(
    mesh,
    thermo
);

// 2. Correct phase fraction, state machine, and thermophysical fields
phaseChange->correct();

// 3. Enthalpy / Energy Equation (solveSolid.H)
fvScalarMatrix hEqn
(
    fvm::ddt(rho, h)
  - fvm::laplacian(thermo.alpha(), h)
 ==
    phaseChange->latentHeatSource()
  + fvm::Sp(phaseChange->latentHeatSp(), h)
  - phaseChange->latentHeatSp() * h_k
);
hEqn.solve();

// 4. Momentum Damping Source (for Enthalpy-Porosity fluid mushy zone flows)
if (phaseChange->active())
{
    UEqn -= phaseChange->momentumSource();
}

// 5. Update state history at completed time-step boundaries
phaseChange->updateHistory();
```

---

## Physics Caveats & Numerical Considerations

> [!IMPORTANT]
> **1. Single-Region Mass Conservation Scope (`massSource() == 0`)**  
> Volumetric Lee models perform in-cell species conversion (conserving cell total density $\rho$). Inter-region mass exchange across coupled boundary interfaces requires conjugate boundary flux coupling.

> [!WARNING]
> **2. Boiling vs. Sub-Boiling Evaporation**  
> $T_{sat}(p)$ evaluates against **total cell pressure** $p$. Phase change occurs when cell temperature exceeds $T_{sat}(p)$. Sub-boiling evaporation driven by partial vapor pressure $p_{v,sat}(T)$ requires species diffusion boundary modeling rather than total-pressure Lee kinetics.

> [!CAUTION]
> **3. Specie Formation Enthalpy ($H_f$) Double Counting**  
> The energy source term $\text{energySource} = - \dot{m} \cdot L$ assumes sensible enthalpy thermo (`sensibleEnthalpy`). If specie thermodynamic properties encode phase latent heat within formation enthalpy $H_f$ (where $H_{f,v} - H_{f,l} = L$), adding explicit $- \dot{m} \cdot L$ will double-count latent heat.

> [!NOTE]
> **4. First-Order $\mathcal{O}(\Delta t)$ Boiling Plateau Temperature Offset**  
> Under large $C_{evap} \Delta t$ and constant volumetric heating $Q$, numerical thermal capping clamps the evaporation rate to $\dot{m} = \frac{\rho C_p \delta}{L \Delta t}$. At steady state ($\dot{m} L = Q$), the plateau temperature lands at $T_{plateau} = T_{sat} + \frac{Q \Delta t}{\rho C_p}$. Halving $\Delta t$ systematically halves this numerical offset.

---

## Utilities, Tutorials, & Verification Suite

### Generating Configuration Files with `writePCDict`

Generate a pre-populated template `phaseChangeDict` with defaults and comments:

```bash
writePCDict
```

### Tutorial Cases

Located in [`tutorials/`](file:///home/lavender/OpenFoamUbu/solvers/phaseChangeFoam/tutorials):

* **[`waterHeaterBoiling`](file:///home/lavender/OpenFoamUbu/solvers/phaseChangeFoam/tutorials/waterHeaterBoiling)**: Multi-region boiling & conjugate heat transfer with fluid Lee phase change.
* **[`pcm_ehc_savE_OM37`](file:///home/lavender/OpenFoamUbu/solvers/phaseChangeFoam/tutorials/pcm_ehc_savE_OM37)**: Thermal energy storage PCM melting/freezing using `EHC`.
* **[`pcm_ep_savE_HS36`](file:///home/lavender/OpenFoamUbu/solvers/phaseChangeFoam/tutorials/pcm_ep_savE_HS36)**: Solidification/melting with fluid flow and Darcy momentum drag using `enthalpyPorosity`.

### Automated Verification Suite

Run the quantitative Python verification suite covering 10 test cases:

```bash
python3 verification/leeModel/runTest.py
```

Covered test cases include analytical ODE rate verification ($< 10^{-7}$ error), Clausius-Clapeyron pressure dependence at 200 kPa, explicit mass and thermal limits, outer-corrector invariance across $nOuterCorrectors \in [1..4]$, heated boiling plateau energy conservation ($0.0000\%$ error), and time-step convergence scaling.
