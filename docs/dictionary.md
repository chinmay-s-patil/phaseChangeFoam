# Dictionary Reference: `phaseChangeDict`

Configuration keywords, subdictionaries, physical units, and numerical controls for `phaseChangeFoam`.

---

## File Placement

In multi-region OpenFOAM setups, each region contains its own directory tree. The phase change dictionary is located at:
```
constant/<regionName>/phaseChangeDict
```
where `<regionName>` matches the region declared in `constant/regionProperties` (e.g. `pcm`, `fluid`, `air`).

---

## Minimal Working Templates

### 1. Solid-Liquid PCM with Effective Heat Capacity (`EHCDarcy`)

Suitable for melting in fluid or solid domains. Incorporates path-integral latent heat and Carman-Kozeny velocity damping:

```foam
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "constant/pcm";
    object      phaseChangeDict;
}

active              true;
meltingModel        EHCDarcy;

forward
{
    Tsolidus        302.68;     // Solidus temperature [K]
    Tliquidus       303.18;     // Liquidus temperature [K]
    latentHeat      80160;      // Latent heat of fusion [J/kg]
}

porosity
{
    Cu              1.6e6;      // Darcy morphology constant [kg/(m3 s)]
    q               1e-3;       // Porosity denominator offset [-]
}
```

### 2. Solid-Liquid PCM with Thermal Hysteresis (`EHC`)

Suitable for solid PCM regions exhibiting supercooling (e.g., salt hydrates or organic waxes):

```foam
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "constant/pcm";
    object      phaseChangeDict;
}

active              true;
meltingModel        EHC;

forward
{
    Tsolidus        309.65;     // Melting start [K]
    Tliquidus       310.65;     // Melting complete [K]
    latentHeat      210000;     // Latent heat [J/kg]
}

reverse
{
    Tsolidus        308.65;     // Freezing complete [K]
    Tliquidus       309.65;     // Freezing start [K]
    latentHeat      210000;     // Latent heat [J/kg]
}

hysteresis
{
    enabled         true;
    subcoolingDelta 1.0;        // Temperature offset [K]
}
```

### 3. Solid-Liquid PCM with Enthalpy-Porosity (`enthalpyPorosity`)

Source-based enthalpy-porosity implementation following Voller & Prakash (1987):

```foam
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "constant/pcm";
    object      phaseChangeDict;
}

active              true;
meltingModel        enthalpyPorosity;

forward
{
    Tsolidus        302.68;     // Solidus temperature [K]
    Tliquidus       303.18;     // Liquidus temperature [K]
    latentHeat      80160;      // Latent heat of fusion [J/kg]
}

porosity
{
    Cu              1.6e6;      // Darcy morphology constant [kg/(m3 s)]
    q               1e-3;       // Porosity denominator offset [-]
}

// Newton projection relaxation factor (keep 1.0 for linear melting window)
lambda              1.0;
convergenceTol      1e-4;
```

### 4. Species Evaporation & Boiling (`Lee`)

Applies volumetric evaporation and condensation to multi-component fluid regions:

```foam
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "constant/fluid";
    object      phaseChangeDict;
}

active              true;
speciesModel        Lee;

liquid              liquid;     // Liquid species name in thermophysicalProperties
vapor               vapor;      // Vapor species name in thermophysicalProperties

C_evap              0.1;        // Evaporation relaxation coefficient [1/s]
C_cond              0.1;        // Condensation relaxation coefficient [1/s]

latentHeat          2.26e6;     // Latent heat of vaporization [J/kg]
Tsat                373.15;     // Saturation temperature at reference pressure [K]

enableTsatP         false;      // Set true to enable Clausius-Clapeyron Tsat(p)
pRef                1e5;        // Reference pressure [Pa]
```

### 5. Fixed Source Verification Baseline (`constantSource`)

```foam
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "constant/fluid";
    object      phaseChangeDict;
}

active              true;
speciesModel        constantSource;

massSource          0.01;       // Net mixture mass source rate [kg/(m3 s)]
energySource        22600.0;    // Energy source rate [W/m3]
```

---

## Keyword Reference Tables

### Top-Level Model Selection
| Keyword | Type | Default | Options | Description |
| :--- | :--- | :--- | :--- | :--- |
| `active` | bool | `true` | `true`, `false` | Enables or disables phase change for this region |
| `meltingModel` | word | none | `EHCDarcy`, `EHC`, `enthalpyPorosity`, `none` | Solid/PCM phase change model |
| `speciesModel` | word | none | `Lee`, `constantSource`, `none` | Fluid species phase change model |

---

### PCM Transition Subdictionaries (`forward` and `reverse`)
| Keyword | Type | Default | Unit | Description |
| :--- | :--- | :--- | :--- | :--- |
| `Tsolidus` | scalar | mandatory | K | Solidus temperature (phase fraction $\alpha = 0$) |
| `Tliquidus` | scalar | mandatory | K | Liquidus temperature (phase fraction $\alpha = 1$) |
| `latentHeat` | scalar | mandatory | J/kg | Latent heat of fusion $L$ |

---

### Hysteresis Subdictionary (`hysteresis`)
| Keyword | Type | Default | Unit | Description |
| :--- | :--- | :--- | :--- | :--- |
| `enabled` | bool | `false` | - | Enables stateful directional hysteresis tracking |
| `subcoolingDelta` | scalar | `0.0` | K | Subcooling offset $T_{\text{melt}} - T_{\text{freeze}}$ |

---

### Porosity Subdictionary (`porosity`)
Active for `EHCDarcy` and `enthalpyPorosity` when operating in fluid regions:
| Keyword | Type | Default | Unit | Description |
| :--- | :--- | :--- | :--- | :--- |
| `Cu` | scalar | `1.6e6` | kg/(m3 s) | Carman-Kozeny morphology coefficient |
| `q` | scalar | `1e-3` | - | Small offset preventing division by zero when $\alpha \to 0$ |

---

### Optional Thermophysical Overrides (`thermophysical`)
Used when distinct solid and liquid thermal conductivities or heat capacities are specified in solid regions:
| Keyword | Type | Default | Unit | Description |
| :--- | :--- | :--- | :--- | :--- |
| `thermoMode` | word | `"none"` | `"none"`, `"custom"` | Enables custom $C_p$ and $k$ properties |
| `Cps` | scalar | none | J/(kg K) | Solid-phase specific heat capacity |
| `Cpl` | scalar | none | J/(kg K) | Liquid-phase specific heat capacity |
| `ks` | scalar | none | W/(m K) | Solid-phase thermal conductivity |
| `kl` | scalar | none | W/(m K) | Liquid-phase thermal conductivity |

---

### Enthalpy-Porosity Specific Settings
| Keyword | Type | Default | Unit | Description |
| :--- | :--- | :--- | :--- | :--- |
| `lambda` | scalar | `1.0` | - | Newton projection relaxation factor |
| `convergenceTol` | scalar | `1e-4` | - | Liquid fraction convergence tolerance |

---

### Fluid Phase Change Settings (`Lee`)
| Keyword | Type | Default | Unit | Description |
| :--- | :--- | :--- | :--- | :--- |
| `liquid` | word | `"liquid"` | - | Name of liquid species in thermophysical database |
| `vapor` | word | `"vapor"` | - | Name of vapor species in thermophysical database |
| `C_evap` | scalar | mandatory | 1/s | Evaporation rate constant |
| `C_cond` | scalar | mandatory | 1/s | Condensation rate constant |
| `latentHeat` | scalar | mandatory | J/kg | Latent heat of vaporization |
| `Tsat` | scalar | `373.15` | K | Reference saturation temperature |
| `enableTsatP` | bool | `false` | - | Enables Clausius-Clapeyron pressure dependence |
| `pRef` | scalar | `1e5` | Pa | Reference saturation pressure |

---

## Solver Numerical Controls

### Multi-Region Coupling & Picard Loop (`system/fvSolution`)

```foam
PIMPLE
{
    nOuterCorrectors     8;          // Multi-region outer iterations (requires >= 2 with phase change)
    outerTolerance       1e-4;       // Temperature convergence tolerance across regions [K]
    nNonLinearCorrectors 25;         // Solid region Picard non-linear corrector passes
    nonLinearTolerance   1e-4;       // Solid region temperature convergence tolerance [K]
    nonLinearToleranceAlpha 1e-4;    // Solid region phase fraction convergence tolerance [-]
}
```

* **Guardrail constraint**: When phase change is active in a fluid region or multi-region domain, `createCoupledRegions.H` enforces `nOuterCorrectors >= 2`. Setting `nOuterCorrectors 1` triggers a `FatalError`.

### Time-Stepping Controls (`system/controlDict`)

```foam
deltaT          0.1;
maxCo           0.5;
maxAlphaCo      0.2;     // Limits time step based on rate of phase fraction change
```

* `maxAlphaCo`: Bounds time step such that cell-wise phase fraction cannot change by more than `maxAlphaCo` within a single step: $\Delta t \le \text{maxAlphaCo} / \max(d\alpha/dt)$, preventing the interface from skipping cells.

---

## Automatic Dictionary Generation (`writePCDict`)

The [`writePCDict`](../src/utilities/writePCDict/writePCDict.C) utility automatically writes formatted `phaseChangeDict` files:

```bash
# Generate interactive EHC dictionary for region 'pcm'
writePCDict -region pcm

# Generate non-interactive Lee model dictionary
writePCDict -region fluid -model Lee -Tsat 373.15 -L 2.26e6 -C_evap 0.1 -C_cond 0.1
```
