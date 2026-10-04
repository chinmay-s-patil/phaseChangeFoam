/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Website:  www.openfoam.com
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

#ifndef phaseChangeModel_H
#define phaseChangeModel_H

#include "fvCFD.H"
#include "solidThermo.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class phaseChangeModel Declaration

 Description
    Unified base class for general multi-region phase-change thermal models
    (solidification/melting, alloys, phase-change materials, ice-water
    transition, and thermal storage) with exact Newton linearization,
    stateful hysteresis, path-integral thermophysics, and fluid momentum source.
\*---------------------------------------------------------------------------*/

class phaseChangeModel
{
protected:

    // Protected data
    const fvMesh& mesh_;
    const solidThermo& thermo_;

    bool active_;

    // Parameters - Melting
    scalar Tlm_; // Melting lower bound [K]
    scalar Tum_; // Melting upper bound [K]
    scalar Lm_;  // Melting latent heat [J/kg]

    // Parameters - Freezing
    scalar Tlf_; // Freezing lower bound [K]
    scalar Tuf_; // Freezing upper bound [K]
    scalar Lf_;  // Freezing latent heat [J/kg]

    // Parameters - Hysteresis
    bool hysteresisActive_;
    scalar reversalTol_; // Temperature hysteresis reversal deadband [K]

    // Parameters - Density
    word densityModel_;
    scalar rhoRef_;
    scalar rhoSolid_;
    scalar rhoLiquid_;

    // Parameters - Thermophysical Mode
    word thermoMode_;
    scalar Cps_;
    scalar Cpl_;
    scalar ks_;
    scalar kl_;

    // Fields
    volScalarField phaseFraction_; // Transformed phase fraction alpha in [0, 1]
    volScalarField phaseState_;    // 0 = solid, 1 = melting, 2 = liquid, 3 = freezing
    volScalarField Cp_;            // Sensible / path-integrated heat capacity [J/(kg K)]
    volScalarField rho_;           // Effective density [kg/m3]
    volScalarField k_;             // Effective thermal conductivity [W/(m K)]
    volScalarField Su_;            // Explicit latent heat source [W/m3]
    volScalarField Sp_;            // Implicit latent heat linearisation coeff [kg/(m3 s)]

    // Trajectory tracking (1.0 = heating, 0.0 = cooling)
    volScalarField heatingTrajectory_;

    // Start-of-step stored fields for committed history
    volScalarField phaseFraction_old_;
    volScalarField heatingTrajectory_old_;

    //- True if phaseFraction was found on disk (or restored from legacy liquidFraction)
    bool phaseFractionRestored_;

    //- Convection suppression flag
    bool suppressConvection_;

    // Protected helper functions
    void readDict();
    void readConvectionDict(const dictionary& pcDict);

    // Exact piecewise path integral helpers
    scalar integralAlphaMelt(scalar T, scalar alpha_old) const;
    scalar integralAlphaFreeze(scalar T, scalar alpha_old) const;

public:

    //- Runtime type information
    TypeName("phaseChangeModel");

    //- Declare runtime construction table
    declareRunTimeSelectionTable
    (
        autoPtr,
        phaseChangeModel,
        dictionary,
        (
            const fvMesh& mesh,
            const solidThermo& thermo
        ),
        (mesh, thermo)
    );

    // Constructors
    phaseChangeModel
    (
        const fvMesh& mesh,
        const solidThermo& thermo
    );

    // Selector
    static autoPtr<phaseChangeModel> New
    (
        const fvMesh& mesh,
        const solidThermo& thermo
    );

    // Destructor
    virtual ~phaseChangeModel() = default;


    // Member Functions

    //- Check if phase change is active for this region
    virtual bool active() const
    {
        return active_;
    }

    //- Correct phase fraction, state machine, thermophysical fields, and Su/Sp
    virtual void correct();

    //- Update trajectory history at completed timestep boundaries
    virtual void updateHistory();

    //- Whether this model suppresses fluid convection
    virtual bool suppressConvection() const
    {
        return suppressConvection_;
    }

    //- Access phase fraction field
    virtual const volScalarField& phaseFraction() const
    {
        return phaseFraction_;
    }

    //- Backward compatibility alias for liquid fraction
    virtual const volScalarField& liquidFraction() const
    {
        return phaseFraction_;
    }

    //- Access phase state field
    virtual const volScalarField& phaseState() const
    {
        return phaseState_;
    }

    //- Access effective heat capacity field
    virtual tmp<volScalarField> CpEff() const
    {
        return Cp_;
    }

    //- Access density field
    virtual tmp<volScalarField> rho() const
    {
        return rho_;
    }

    //- Access thermal conductivity field
    virtual tmp<volScalarField> k() const
    {
        return k_;
    }

    //- Access explicit latent heat source term
    virtual tmp<volScalarField> latentHeatSource() const
    {
        return Su_;
    }

    //- Access implicit latent heat linearisation coefficient
    virtual tmp<volScalarField> latentHeatSp() const
    {
        return Sp_;
    }

    //- Access momentum source term (default: zero, override in models with Darcy damping)
    virtual tmp<volVectorField> momentumSource() const
    {
        return tmp<volVectorField>::New
        (
            IOobject("momentumSource", mesh_.time().timeName(), mesh_),
            mesh_,
            dimensionedVector("zero", dimForce/dimVolume, Zero)
        );
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
