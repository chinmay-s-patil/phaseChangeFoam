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

To implement a new solid-liquid phase change formulation (e.g. `myPCMModel`):

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
    virtual tmp<volScalarField> kEff() const override;

    //- Explicit latent heat source term [W/m3]
    virtual tmp<volScalarField> Su() const override;

    //- Implicit linearisation coefficient [kg/(m3 s)]
    virtual tmp<volScalarField> Sp() const override;

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
    myParam_(subOrEmptyDict("myPCMCoeffs").getOrDefault<scalar>("myParam", 1.0))
{}

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

tmp<volScalarField> myPCMModel::kEff() const
{
    return phaseChangeModel::kEff();
}

tmp<volScalarField> myPCMModel::Su() const
{
    return phaseChangeModel::Su();
}

tmp<volScalarField> myPCMModel::Sp() const
{
    return phaseChangeModel::Sp();
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

    //- Net mixture energy source rate per unit volume [W/m3]
    virtual tmp<volScalarField> energySource() const override;

    //- Species mass production rate [kg/(m3 s)]
    virtual tmp<volScalarField> speciesSource(const word& speciesName) const override;

    //- Update rate fields and kinetic state
    virtual void correct() override;

    //- Diagnostics check for mass conservation
    virtual void checkMassConservation() const override;
};

} // End namespace Foam

#endif
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
