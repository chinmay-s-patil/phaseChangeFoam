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

All phase change configurations require a top-level `phaseChange` subdictionary in `constant/<regionName>/phaseChangeDict`. The solver validates keys using `checkAllowedKeys()`.

### 1. Solid-Liquid PCM with Effective Heat Capacity (`EHCDarcy`)

Suitable for melting in fluid or solid domains. Incorporates path-based sensible heat capacity, linearised latent heat sources, and Carman-Kozeny velocity damping:

```foam
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "constant/pcm";
    object      phaseChangeDict;
}

active          true;

phaseChange
{
    active          true;
    phaseChangeMode EHCDarcy;

    forward
    {
        T_lowerBound    302.68;     // Solidus temperature [K]
        T_upperBound    303.18;     // Liquidus temperature [K]
        latentHeat      80160;      // Latent heat of fusion [J/kg]
    }

    porosity
    {
        Cu              1.6e6;      // Benchmark value (default: 1e5) [kg/(m3 s)]
        q               1e-3;       // Porosity denominator offset [-]
    }
}
```

### 2. Solid-Liquid PCM with Thermal Hysteresis (`EHC`)

Suitable for solid PCM regions or pure-conduction domains exhibiting supercooling (e.g., salt hydrates or organic waxes). Does not include Darcy momentum drag; for fluid regions with convection, use `EHCDarcy` (or set `convection { suppress true; }` to run pure conduction in a fluid domain):

```foam
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "constant/pcm";
    object      phaseChangeDict;
}

active          true;

phaseChange
{
    active          true;
    phaseChangeMode EHC;
    direction       both;

    forward
    {
        T_lowerBound    309.65;     // Melting start [K]
        T_upperBound    310.65;     // Melting complete [K]
        latentHeat      210000;     // Latent heat of fusion [J/kg]
    }

    reverse
    {
        T_lowerBound    308.65;     // Freezing complete [K]
        T_upperBound    309.65;     // Freezing start [K]
        latentHeat      210000;     // Latent heat of freezing [J/kg]
    }

    hysteresis
    {
        active              true;   // Enable path-dependent thermal hysteresis
        reversalTolerance   1e-3;   // Reversal detection temperature tolerance [K]
    }
}
```

### 3. Solid-Liquid PCM with Enthalpy-Porosity (`enthalpyPorosity`)

Source-based enthalpy-porosity implementation following Voller & Prakash (1987). Supports a single phase-transition interval (`forward` melting window) and does not support separate melting/freezing curves, reverse transitions, or stateful thermal hysteresis (configuring `reverse` or `hysteresis` subdictionaries, or setting `direction` other than `both`, triggers a runtime fatal error):

```foam
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "constant/pcm";
    object      phaseChangeDict;
}

active          true;

phaseChange
{
    active          true;
    phaseChangeMode enthalpyPorosity;

    forward
    {
        T_lowerBound    302.68;     // Solidus temperature [K]
        T_upperBound    303.18;     // Liquidus temperature [K]
        latentHeat      80160;      // Latent heat of fusion [J/kg]
    }

    porosity
    {
        Cu              1.6e6;      // Benchmark value (default: 1e5) [kg/(m3 s)]
        q               1e-3;       // Porosity denominator offset [-]
    }

    lambda          1.0;        // Under-relaxation factor in (0, 1]
    convergenceTol  1e-4;       // Convergence tolerance on max|Delta f|
}
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

active          true;

phaseChange
{
    active          true;
    type            Lee;

    liquid          liquid;     // Liquid species name in thermophysicalProperties
    vapor           vapor;      // Vapor species name in thermophysicalProperties

    C_evap          0.1;        // Evaporation relaxation coefficient [1/s]
    C_cond          0.1;        // Condensation relaxation coefficient [1/s]

    latentHeat      2.26e6;     // Latent heat of vaporization [J/kg]
    Tsat            373.15;     // Saturation temperature at reference pressure [K]
    pRef            101325;     // Reference pressure [Pa]
    enableTsatP     false;      // Set true to enable Clausius-Clapeyron Tsat(p)
}
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

active          true;

phaseChange
{
    active          true;
    type            constantSource;

    speciesName     H2O;        // Name of species receiving mass source
    massSource      0.01;       // Net mixture mass source rate [kg/(m3 s)]
    energySource    22600.0;    // Energy source rate [W/m3]
}
```

---

## Keyword Reference Tables

### Top-Level and `phaseChange` Model Selection
| Keyword | Scope | Type | Default | Options | Description |
| :--- | :--- | :--- | :--- | :--- | :--- |
| `active` | Root & `phaseChange` | bool | `true` | `true`, `false` | Enables or disables phase change for this region |
| `phaseChangeMode` / `meltingModel` | `phaseChange` | word | none | `EHCDarcy`, `EHC`, `enthalpyPorosity`, `none` | Solid/PCM phase change model |
| `type` / `speciesModel` | `phaseChange` | word | none | `Lee`, `constantSource`, `none` | Fluid species phase change model |
| `direction` | `phaseChange` | word | `both` | `both`, `forward`, `reverse` | Transition directionality for solid PCM (`enthalpyPorosity` supports only `both`) |

---

### PCM Transition Subdictionaries (`forward` and `reverse`)
The `forward` block specifies the primary phase change window (melting). The `reverse` block is optional and supported by `EHC`/`EHCDarcy` for separate freezing paths (`enthalpyPorosity` supports only `forward`).
| Keyword | Type | Default | Unit | Description |
| :--- | :--- | :--- | :--- | :--- |
| `T_lowerBound` | scalar | mandatory | K | Solidus / lower bound temperature (phase fraction $\alpha = 0$) |
| `T_upperBound` | scalar | mandatory | K | Liquidus / upper bound temperature (phase fraction $\alpha = 1$) |
| `latentHeat` | scalar | mandatory | J/kg | Latent heat of fusion / freezing $L$ |

---

### Hysteresis Subdictionary (`hysteresis`)
Supported by `EHC` and `EHCDarcy` (rejected by `enthalpyPorosity`):
| Keyword | Type | Default | Unit | Description |
| :--- | :--- | :--- | :--- | :--- |
| `active` | bool | `false` | - | Enables stateful directional hysteresis tracking |
| `reversalTolerance` | scalar | `1e-3` | K | Temperature reversal detection threshold |

---

### Porosity Subdictionary (`porosity`)
Active for `EHCDarcy` and `enthalpyPorosity` when operating in fluid regions:
| Keyword | Type | Default | Unit | Description |
| :--- | :--- | :--- | :--- | :--- |
| `Cu` (or `A_cu`) | scalar | `1e5` | kg/(m3 s) | Carman-Kozeny morphology coefficient / Darcy resistance constant (e.g. `1.6e6` in benchmark configurations) |
| `q` (or `eps`) | scalar | `1e-3` | - | Small offset preventing division by zero when $\alpha \to 0$ |

---

### Convection Suppression Subdictionary (`convection`)
Optional subdictionary to control fluid velocity coupling and Darcy momentum drag:
| Keyword | Type | Default | Unit | Description |
| :--- | :--- | :--- | :--- | :--- |
| `suppress` | bool | `true` (solid/no $\mathbf{U}$), `false` (`EHCDarcy`/`enthalpyPorosity`) | - | Suppresses fluid velocity coupling and Darcy momentum drag in the energy and momentum equations. In `EHC`, requires `true` if instantiated in a fluid region with velocity field $\mathbf{U}$. In pure solid/conduction regions without $\mathbf{U}$, defaults to `true` (setting `false` throws a `FatalIOError`). In `EHCDarcy`, setting `suppress true` freezes flow coupling while retaining the model formulation. |

---

### Optional Thermophysical Overrides (`thermophysical`)
Used when distinct solid and liquid thermal conductivities or heat capacities are specified in solid regions:
| Keyword | Type | Default | Unit | Description |
| :--- | :--- | :--- | :--- | :--- |
| `mode` | word | `"thermo"` | `"thermo"`, `"custom"` | Enables custom $C_p$ and $k$ properties |
| `CpSolid` (or `Cps`) | scalar | none | J/(kg K) | Solid-phase specific heat capacity |
| `CpLiquid` (or `Cpl`) | scalar | none | J/(kg K) | Liquid-phase specific heat capacity |
| `kSolid` (or `ks`) | scalar | none | W/(m K) | Solid-phase thermal conductivity |
| `kLiquid` (or `kl`) | scalar | none | W/(m K) | Liquid-phase thermal conductivity |

---

### Optional Density Overrides (`density`)
| Keyword | Type | Default | Unit | Description |
| :--- | :--- | :--- | :--- | :--- |
| `model` | word | `"thermo"` | `"thermo"`, `"linear"` | Density formulation (`thermo`: constant conservative density) |
| `rhoRef` | scalar | `1000.0` | kg/m3 | Reference density |
| `rhoSolid` | scalar | `rhoRef` | kg/m3 | Solid-phase density |
| `rhoLiquid` | scalar | `rhoRef` | kg/m3 | Liquid-phase density |
| `allowNonConservativeDensity` | bool | `false` | - | Required `true` if $\rho_{\text{solid}} \ne \rho_{\text{liquid}}$ in solid regions |

---

### Enthalpy-Porosity Specific Settings
| Keyword | Type | Default | Unit | Description |
| :--- | :--- | :--- | :--- | :--- |
| `lambda` | scalar | `1.0` | - | Newton projection relaxation factor in range $(0, 1]$ |
| `convergenceTol` | scalar | `0.01` | - | Liquid fraction convergence tolerance |

---

### Fluid Phase Change Settings (`Lee`)
| Keyword | Type | Default | Unit | Description |
| :--- | :--- | :--- | :--- | :--- |
| `liquid` | word | `"liquid"` | - | Name of liquid species in thermophysical database |
| `vapor` | word | `"vapor"` | - | Name of vapor species in thermophysical database |
| `C_evap` | scalar | `0.1` | 1/s | Evaporation rate constant |
| `C_cond` | scalar | `0.1` | 1/s | Condensation rate constant |
| `latentHeat` | scalar | mandatory | J/kg | Latent heat of vaporization |
| `Tsat` (or `T_sat`) | scalar | `373.15` | K | Reference saturation temperature |
| `enableTsatP` | bool | `false` | - | Enables Clausius-Clapeyron pressure dependence |
| `pRef` | scalar | `101325` | Pa | Reference saturation pressure |

---

### Fluid Phase Change Settings (`constantSource`)
| Keyword | Type | Default | Unit | Description |
| :--- | :--- | :--- | :--- | :--- |
| `speciesName` | word | `"H2O"` | - | Name of species receiving mass source |
| `massSource` | scalar | `0.0` | kg/(m3 s) | Net mixture mass source rate (or `speciesSource`) |
| `energySource` | scalar | `-1e6` | W/m3 | Energy source rate |

---

## Solver Numerical Controls

### Multi-Region Coupling & Picard Loop (`system/fvSolution`)

```foam
coupledEnergyField   false;      // Enable monolithic implicit energy matrix coupling across interfaces

PIMPLE
{
    nOuterCorrectors     8;          // Multi-region outer iterations (requires >= 2 with phase change)
    outerTolerance       1e-4;       // Temperature convergence tolerance across regions [K]
    nNonLinearCorrectors 25;         // Solid region Picard non-linear corrector passes
    nonLinearTolerance   1e-4;       // Solid region temperature convergence tolerance [K]
    nonLinearToleranceAlpha 1e-4;    // Solid region phase fraction convergence tolerance [-]
    frozenFlow           false;      // Freeze velocity/momentum solve and only advance energy & phase change
}
```

* **Guardrail constraint**: When phase change is active in a fluid region or multi-region domain, `createCoupledRegions.H` enforces `nOuterCorrectors >= 2`. Setting `nOuterCorrectors 1` triggers a `FatalError`.
* **`coupledEnergyField`**: Optional boolean setting at the root of `system/fvSolution` (default: `false`). When set to `true`, or when any conjugate patch sets `useImplicit true;`, `phaseChangeMultiRegionFoam` activates monolithic implicit matrix assembly across solid-fluid and solid-solid interfaces (`createCoupledRegions.H`). The energy equations of the coupled regions are assembled into a unified block matrix (`fvMatrixAssemblyPtr`) and solved simultaneously, eliminating Picard interface iteration lag across conjugate heat transfer boundaries.
* **`frozenFlow`**: Optional boolean setting in fluid `PIMPLE` or `SIMPLE` subdictionaries (default: `false`). When set to `true`, fluid momentum, pressure correction, continuity, and turbulence equations (`UEqn`, `pEqn`, `rhoEqn`) are bypassed. The solver holds velocity $\mathbf{U}$ and mass flux $\phi$ static, advancing only thermal energy (`EEqn`), species transport (`YEqn`), and phase change.

### Time-Stepping Controls (`system/controlDict`)

```foam
deltaT          0.1;
maxCo           0.5;
maxDi           10.0;    // Limits time step based on solid thermal diffusion number
maxAlphaCo      0.2;     // Limits time step based on rate of phase fraction change
```

* **`maxDi`**: Maximum allowable diffusion number in solid/conduction regions ($\text{Di} = \frac{k}{\rho c_p} \frac{\Delta t}{\Delta x^2}$), read from `controlDict` (default: `10.0`). In multi-region transient runs with adjustable time stepping, the time step is bounded by both fluid Courant limits and solid diffusion: $\Delta t_{\text{solid}} \le \text{maxDi} / (\text{Di}_{\text{num}} + \epsilon)$.
* **`maxAlphaCo`**: Bounds time step such that cell-wise phase fraction cannot change by more than `maxAlphaCo` within a single step: $\Delta t \le \text{maxAlphaCo} / \max(d\alpha/dt)$, preventing the interface from skipping cells.

---

## Automatic Dictionary Generation (`writePCDict`)

The [`writePCDict`](../src/utilities/writePCDict/writePCDict.C) utility automatically writes formatted `phaseChangeDict` files with all valid subdictionaries and defaults:

```bash
# Generate default template in constant/phaseChange/phaseChangeDict
writePCDict

# Generate template for region 'pcm' in constant/pcm/phaseChangeDict
writePCDict -region pcm

# Generate template for fluid region 'fluid' in constant/fluid/phaseChangeDict
writePCDict -region fluid
```

> [!NOTE]
> `writePCDict` outputs a template populated with baseline defaults across solid and fluid models. Generated dictionaries must be reviewed and populated with the appropriate material properties, operating temperatures, and domain-specific models for your simulation case.
