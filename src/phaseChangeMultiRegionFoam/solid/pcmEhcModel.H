/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Website:  www.openfoam.com
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

#ifndef pcmEhcModel_H
#define pcmEhcModel_H

#include "pcmPhaseChangeModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class pcmEhcModel Declaration
\*---------------------------------------------------------------------------*/

class pcmEhcModel
:
    public pcmPhaseChangeModel
{
private:

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

    // Parameters - Density
    word densityModel_;
    scalar rhoRef_;
    scalar rhoSolid_;
    scalar rhoLiquid_;

    // Parameters - Thermophysical
    word thermoMode_;
    scalar Cps_;
    scalar Cpl_;
    scalar ks_;
    scalar kl_;

    // Calculated Fields
    volScalarField CpEff_;
    volScalarField rho_;
    volScalarField k_;

    // Trajectory tracking (1.0 = heating, 0.0 = cooling)
    volScalarField heatingTrajectory_;

    // Start-of-step stored fields for committed alphaL_old formulation across Picard iterations
    volScalarField liquidFraction_old_;
    volScalarField heatingTrajectory_old_;

    // Helper functions
    void readDict();

    // Exact piecewise path integral helpers for sensible heat capacity integration
    scalar integralAlphaMelt(scalar T, scalar alpha_old) const;
    scalar integralAlphaFreeze(scalar T, scalar alpha_old) const;

public:

    //- Runtime type information
    TypeName("EHC");

    // Constructors
    pcmEhcModel
    (
        const fvMesh& mesh,
        const solidThermo& thermo
    );

    // Destructor
    virtual ~pcmEhcModel() = default;


    // Member Functions

    virtual bool active() const override
    {
        return active_;
    }

    virtual void correct() override;

    virtual void updateHistory() override;

    virtual tmp<volScalarField> CpEff() const override
    {
        return CpEff_;
    }

    virtual tmp<volScalarField> rho() const override
    {
        return rho_;
    }

    virtual tmp<volScalarField> k() const override
    {
        return k_;
    }

    virtual tmp<volScalarField> latentHeatSource() const override
    {
        return tmp<volScalarField>::New
        (
            IOobject("latentHeatSource", mesh_.time().timeName(), mesh_),
            mesh_,
            dimensionedScalar("zero", dimEnergy/dimTime/dimVolume, 0.0)
        );
    }

    virtual tmp<volVectorField> momentumSource() const override
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
