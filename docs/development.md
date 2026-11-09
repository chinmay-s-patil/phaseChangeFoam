# Developer & Model Integration Guide

Architecture, execution sequence, numerical coupling, and step-by-step instructions for extending `phaseChangeFoam`.

---

## Source Tree Structure

```
src/
├── solvers/
│   ├── phaseChangeMultiRegionFoam/        # Transient multi-region CHT & phase change
│   │   ├── phaseChangeMultiRegionFoam.C   # Main application driver & PIMPLE loop
│   │   ├── fluid/                         # Fluid region models and equation solvers
│   │   │   ├── fluidPhaseChangeModel.[CH] # Base class for fluid species phase change
│   │   │   ├── leeFluidPhaseChangeModel.[CH]
│   │   │   ├── constantSourceFluidPhaseChangeModel.[CH]
│   │   │   ├── noneFluidPhaseChangeModel.[CH]
│   │   │   ├── solveFluid.H               # Fluid corrector sequence
│   │   │   ├── UEqn.H                     # Momentum with implicit Darcy drag
│   │   │   ├── EEqn.H                     # Sensible enthalpy equation
│   │   │   ├── pEqn.H                     # Pressure-continuity equation
│   │   │   └── YEqn.H                     # Species mass fractions
│   │   ├── solid/                         # Solid & PCM models and equation solvers
│   │   │   ├── phaseChangeModel.[CH]      # Base class for solid-liquid PCM
│   │   │   ├── ehcPhaseChangeModel.[CH]
│   │   │   ├── ehcDarcyPhaseChangeModel.[CH]
│   │   │   ├── enthalpyPorosityPhaseChangeModel.[CH]
│   │   │   ├── nonePhaseChangeModel.[CH]
│   │   │   └── solveSolid.H               # Non-linear Picard loop
│   │   └── include/                       # Coupling and time-stepping headers
│   │       ├── createCoupledRegions.H     # fvMatrixAssembly initialization
│   │       ├── correctThermos.H           # Post-solve thermo update
│   │       └── setMultiRegionDeltaT.H     # Adaptive time step with maxAlphaCo
│   └── phaseChangeMultiRegionSimpleFoam/  # Steady-state multi-region CHT solver
└── utilities/
    └── writePCDict/                       # Dictionary generation utility
```

---

## Compilation

Build the entire solver suite using the top-level script:
```bash
./Allwmake
```
To remove compiled object files and binaries:
```bash
./Allwclean
```

---

## Solver Execution Sequence

### PIMPLE Outer Loop Lifecycle (`phaseChangeMultiRegionFoam.C`)

1. **Adaptive Time Step Evaluation**:
   * Evaluates Courant numbers across fluid regions (`compressibleMultiRegionCourantNo.H`).
   * Evaluates diffusion numbers across solid regions (`solidRegionDiffusionNo.H`).
   * Limits $\Delta t$ via `setMultiRegionDeltaT.H` ensuring both Courant bounding and phase-fraction Courant bounding:
     $$\Delta t \le \frac{\text{maxAlphaCo}}{\max(d\alpha/dt)}$$

2. **PIMPLE Outer Loop** (`for (int oCorr=0; oCorr<nOuterCorr; ++oCorr)`):
   * **Fluid regions (`solveFluid.H`)**:
     * If $oCorr == 0$: update species evaporation/condensation rates (`phaseChange.correct()`).
     * If PCM model is active in fluid: update phase fraction and apparent capacity (`pcmModel.correct()`).
     * Assemble and solve momentum (`UEqn.H`), incorporating implicit Carman-Kozeny Darcy drag $\mathbf{S}_u$.
     * Assemble and solve species mass fractions (`YEqn.H`).
     * Assemble sensible enthalpy equation (`EEqn.H`).
     * If uncoupled: solve PISO pressure-velocity correctors (`pEqn.H`).
   * **Solid regions (`solveSolid.H`)**:
     * Executes non-linear Picard iteration loop (`nNonLinearCorrectors`) to resolve the non-linear coupling between temperature $T$ and apparent heat capacity $C_{p,\text{eff}}(T)$ or explicit source $S_h$.
   * **Coupled energy matrix solution**:
     * When `coupled` is active: `fvMatrixAssemblyPtr->solve()` solves the monolithic global energy matrix across all solid and fluid regions simultaneously.
     * Evaluates outer temperature change: $\max(\Delta T_{\text{outer}})$.
     * **Early convergence handling**: If $\max(\Delta T_{\text{outer}}) < \text{outerTolerance}$, outer iterations are declared converged. Before breaking out of the outer loop, the solver executes a final PISO pressure-velocity corrector loop across all fluid regions (`pEqn.H`), ensuring the velocity and pressure fields are fully consistent with the converged temperature field.
     * If not early converged: executes regular PISO loop (`pEqn.H`) across fluid regions.

---

## Implementing a New Phase Change Model

Models are instantiated via OpenFOAM runtime selection tables (`runTimeSelectionTable`).

### 1. Adding a Solid / PCM Model (`phaseChangeModel`)

> [!IMPORTANT]
> **Custom PCM Model Architecture & Implementation Limitation**:
> The current `phaseChangeModel` base class performs model-name validation and dictionary-key validation internally. Registering a new class in the runtime selection table (`addToRunTimeSelectionTable`) is therefore not sufficient by itself to introduce an independent PCM formulation:
> 1. **Model Name Validation**: In `phaseChangeModel::readDict()`, the solver checks the configured `phaseChangeMode` against `"ehc"`, `"ehcdarcy"`, and `"enthalpyporosity"`. An unrecognized model name will be marked inactive (`active_ = false`) by the base class unless the base-class selection and dictionary-reading logic in [`src/solvers/phaseChangeMultiRegionFoam/solid/phaseChangeModel.C`](../src/solvers/phaseChangeMultiRegionFoam/solid/phaseChangeModel.C) is extended to recognise the new model.
> 2. **Dictionary Key Validation**: `phaseChangeModel::readDict()` validates keys in `phaseChangeDict` via `checkAllowedKeys()`. Custom subdictionaries (such as `myPCMCoeffs`) must be registered in the allowed-key list.
> 3. **Base Implementation Reuse vs. Independent Formulations**: In the code example below, delegation calls such as `phaseChangeModel::correct()`, `phaseChangeModel::CpEff()`, and `phaseChangeModel::latentHeatSource()` reuse the existing base implementation (which computes path-secant heat capacity and thermal hysteresis). A genuinely different formulation must implement its own phase-change state updates and provide appropriate field contributions rather than delegating to the base class.
>
> The following example illustrates the class structure and available virtual extension points; it is an architectural reference rather than a drop-in implementation.

#### Header: `myPCMModel.H`
```cpp
#ifndef myPCMModel_H
#define myPCMModel_H

#include "phaseChangeModel.H"

namespace Foam
{

class myPCMModel
:
    public phaseChangeModel
{
    // Model-specific parameters
    scalar myParam_;

public:

    TypeName("myPCM");

    // Constructor
    myPCMModel(const fvMesh& mesh, const basicThermo& thermo);

    // Destructor
    virtual ~myPCMModel() = default;

    //- Update phase fraction and state variables
    virtual void correct() override;

    //- Update committed history at completed time steps
    virtual void updateHistory() override;

    //- Hook called after energy solve (for explicit projection methods)
    virtual void postEnergySolve(volScalarField& he, const bool finalIter = false) override;

    //- Apply effective conductivity to thermo database
    virtual void applyKappaToThermo() override;

    //- Effective heat capacity field [J/(kg K)]
    virtual tmp<volScalarField> CpEff() const override;

    //- Effective density field [kg/m3]
    virtual tmp<volScalarField> rhoEff() const override;

    //- Effective thermal conductivity field [W/(m K)]
    virtual tmp<volScalarField> k() const override;

    //- Explicit latent heat source term [W/m3]
    virtual tmp<volScalarField> latentHeatSource() const override;

    //- Implicit linearisation coefficient [kg/(m3 s)]
    virtual tmp<volScalarField> latentHeatSp() const override;

    //- Darcy momentum drag coefficient [kg/(m3 s)] (for fluid regions)
    virtual tmp<volScalarField> momentumSp() const override;
};

} // End namespace Foam

#endif
```

#### Source: `myPCMModel.C`
```cpp
#include "myPCMModel.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug(myPCMModel, 0);
addToRunTimeSelectionTable(phaseChangeModel, myPCMModel, dictionary);

myPCMModel::myPCMModel(const fvMesh& mesh, const basicThermo& thermo)
:
    phaseChangeModel(mesh, thermo),
    myParam_(1.0)
{
    // Read model-specific parameters from constant/<regionName>/phaseChangeDict
    IOobject dictIO
    (
        "phaseChangeDict",
        mesh_.time().constant(),
        mesh_,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    );

    if (dictIO.typeHeaderOk<dictionary>(true))
    {
        IOdictionary dict(dictIO);
        if (dict.found("phaseChange"))
        {
            const dictionary& pcDict = dict.subDict("phaseChange");
            if (pcDict.found("myPCMCoeffs"))
            {
                const dictionary& modelDict = pcDict.subDict("myPCMCoeffs");
                myParam_ = modelDict.getOrDefault<scalar>("myParam", 1.0);
            }
        }
    }
}

void myPCMModel::correct()
{
    phaseChangeModel::correct();
    // Update model-specific field variables here
}

void myPCMModel::updateHistory()
{
    phaseChangeModel::updateHistory();
}

void myPCMModel::postEnergySolve(volScalarField& he, const bool finalIter)
{}

void myPCMModel::applyKappaToThermo()
{
    phaseChangeModel::applyKappaToThermo();
}

tmp<volScalarField> myPCMModel::CpEff() const
{
    return phaseChangeModel::CpEff();
}

tmp<volScalarField> myPCMModel::rhoEff() const
{
    return phaseChangeModel::rhoEff();
}

tmp<volScalarField> myPCMModel::k() const
{
    return phaseChangeModel::k();
}

tmp<volScalarField> myPCMModel::latentHeatSource() const
{
    return phaseChangeModel::latentHeatSource();
}

tmp<volScalarField> myPCMModel::latentHeatSp() const
{
    return phaseChangeModel::latentHeatSp();
}

tmp<volScalarField> myPCMModel::momentumSp() const
{
    return phaseChangeModel::momentumSp();
}

} // End namespace Foam
```

---

### 2. Adding a Fluid Phase Change Model (`fluidPhaseChangeModel`)

To implement a new fluid-region phase change or evaporation model:

#### Header: `myFluidPhaseChangeModel.H`
```cpp
#ifndef myFluidPhaseChangeModel_H
#define myFluidPhaseChangeModel_H

#include "fluidPhaseChangeModel.H"

namespace Foam
{

class myFluidPhaseChangeModel
:
    public fluidPhaseChangeModel
{
    scalar myRateConst_;

public:

    TypeName("myFluidModel");

    // Constructor
    myFluidPhaseChangeModel
    (
        const fvMesh& mesh,
        const rhoReactionThermo& thermo,
        const volVectorField& U,
        const surfaceScalarField& phi
    );

    virtual ~myFluidPhaseChangeModel() = default;

    //- Net mixture mass source rate per unit volume [kg/(m3 s)]
    virtual tmp<volScalarField> massSource() const override;

    //- Species mass production rate by specie index [kg/(m3 s)] (pure virtual)
    virtual tmp<volScalarField> speciesSource(const label specieIndex) const override;

    //- Species mass production rate by specie name [kg/(m3 s)] (optional overload)
    virtual tmp<volScalarField> speciesSource(const word& specieName) const override;

    //- Net mixture energy source rate per unit volume [W/m3]
    virtual tmp<volScalarField> energySource() const override;

    //- Update rate fields and kinetic state
    virtual void correct() override;

    //- Diagnostics check for mass conservation: returns true if sum(speciesSource(i)) == massSource()
    virtual bool checkMassConservation
    (
        const scalar relTol = 1e-4,
        const scalar absTol = 1e-10
    ) const override;
};

} // End namespace Foam

#endif
```

#### Source: `myFluidPhaseChangeModel.C`
```cpp
#include "myFluidPhaseChangeModel.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug(myFluidPhaseChangeModel, 0);
addToRunTimeSelectionTable(fluidPhaseChangeModel, myFluidPhaseChangeModel, dictionary);

myFluidPhaseChangeModel::myFluidPhaseChangeModel
(
    const fvMesh& mesh,
    const rhoReactionThermo& thermo,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    fluidPhaseChangeModel(mesh, thermo, U, phi),
    myRateConst_(dict_.getOrDefault<scalar>("myRateConst", 1.0))
{}

tmp<volScalarField> myFluidPhaseChangeModel::massSource() const
{
    return tmp<volScalarField>::New
    (
        IOobject("massSource", mesh_.time().timeName(), mesh_),
        mesh_,
        dimensionedScalar("zero", dimMass/dimVolume/dimTime, 0.0)
    );
}

tmp<volScalarField> myFluidPhaseChangeModel::speciesSource(const label specieIndex) const
{
    return tmp<volScalarField>::New
    (
        IOobject("speciesSource", mesh_.time().timeName(), mesh_),
        mesh_,
        dimensionedScalar("zero", dimMass/dimVolume/dimTime, 0.0)
    );
}

tmp<volScalarField> myFluidPhaseChangeModel::speciesSource(const word& specieName) const
{
    return fluidPhaseChangeModel::speciesSource(specieName);
}

tmp<volScalarField> myFluidPhaseChangeModel::energySource() const
{
    return tmp<volScalarField>::New
    (
        IOobject("energySource", mesh_.time().timeName(), mesh_),
        mesh_,
        dimensionedScalar("zero", dimEnergy/dimVolume/dimTime, 0.0)
    );
}

void myFluidPhaseChangeModel::correct()
{
    // Update kinetic rates and mass/energy transfer fields
}

bool myFluidPhaseChangeModel::checkMassConservation
(
    const scalar relTol,
    const scalar absTol
) const
{
    return fluidPhaseChangeModel::checkMassConservation(relTol, absTol);
}

} // End namespace Foam
```

---

## Compiling & Testing New Models

1. Add the new `.C` files to `Make/files`:
   ```
   fluid/myFluidPhaseChangeModel.C
   solid/myPCMModel.C
   ```
2. Recompile:
   ```bash
   ./Allwmake
   ```
3. Run the automated verification test harness:
   ```bash
   python3 verification/runAll.py
   python3 validation/runAll.py --fast
   ```
   Add a targeted test case under `verification/` verifying analytical energy conservation or mass transfer invariance for the new model before deploying to production.
