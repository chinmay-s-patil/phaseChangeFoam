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
public:

    //- Phase change direction mode restriction
    enum class direction
    {
        both,
        forward,
        reverse
    };

    static const Enum<direction> directionNames;

protected:

    // Protected data
    const fvMesh& mesh_;
    const solidThermo& thermo_;

    bool active_;

    // Direction mode restriction
    direction dir_;

    // Parameters - Melting / Forward
    scalar Tlm_; // Melting/Forward lower bound [K]
    scalar Tum_; // Melting/Forward upper bound [K]
    scalar Lm_;  // Melting/Forward latent heat [J/kg]

    // Parameters - Freezing / Reverse
    scalar Tlf_; // Freezing/Reverse lower bound [K]
    scalar Tuf_; // Freezing/Reverse upper bound [K]
    scalar Lf_;  // Freezing/Reverse latent heat [J/kg]

    // Parameters - Hysteresis
    bool hysteresisActive_;
    scalar reversalTol_; // Temperature hysteresis reversal deadband [K]

    // Parameters - Density
    word densityModel_;
    scalar rhoRef_;
    scalar rhoSolid_;
    scalar rhoLiquid_;
    bool allowNonConservativeDensity_;

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
    volScalarField T_reversal_;

    // Start-of-step stored fields for committed history
    volScalarField phaseFraction_old_;
    volScalarField heatingTrajectory_old_;
    volScalarField T_reversal_old_;

    //- True if phaseFraction was found on disk (or restored from legacy liquidFraction)
    bool phaseFractionRestored_;

    //- Convection suppression flag
    bool suppressConvection_;

    // Protected helper structs and functions
    struct StateResult
    {
        scalar traj;
        scalar Trev;
        scalar aVal;
        scalar stateVal;
        bool isPlateau;
        scalar dAlpha_dT;
    };

    StateResult evalState
    (
        scalar T,
        scalar Told,
        scalar aOld,
        scalar trajOld,
        scalar TrevOld
    ) const;

    void readDict();
    void readConvectionDict(const dictionary& pcDict);

    // Exact sensible enthalpy integral for custom thermo mode
    scalar hSens(scalar T, scalar traj) const;

    // Exact step sensible enthalpy change considering hysteresis plateaus
    scalar deltaHSens(scalar Tnew, scalar Told, scalar aOld, scalar trajNew, scalar trajOld) const;

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
        const solidThermo& thermo,
        bool suppressConvection = true
    );
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

    //- Access density field for energy equation (returns background thermo.rho() for reference independence)
    virtual tmp<volScalarField> rho() const
    {
        return thermo_.rho();
    }

    //- Access physical phase density field rho(alpha)
    virtual tmp<volScalarField> rhoEff() const
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
