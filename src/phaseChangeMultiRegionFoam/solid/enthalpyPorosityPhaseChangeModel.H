/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Website:  www.openfoam.com
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

#ifndef pcmEnthalpyPorosityModel_H
#define pcmEnthalpyPorosityModel_H

#include "pcmPhaseChangeModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class pcmEnthalpyPorosityModel Declaration

 Description
    Enthalpy-Porosity PCM model for phase-change simulation.
    Note: Full classical enthalpy-porosity with momentum damping (Darcy source)
    requires a fluid region where UEqn is solved. In solid PCM regions, only
    sensible heat conduction and latent heat source terms are active.
\*---------------------------------------------------------------------------*/

class pcmEnthalpyPorosityModel
:
    public pcmPhaseChangeModel
{
private:

    // Parameters - Melting & Freezing
    scalar Tlm_;
    scalar Tum_;
    scalar Lm_;

    // Parameters - Darcy & Buoyancy
    scalar beta_; // Thermal expansion coefficient [1/K]
    scalar Cu_;   // Darcy penalty coefficient [kg/(m3 s)]
    scalar q_;    // Smoothing constant

    // Parameters - Density & Thermo
    scalar rhoRef_;
    scalar rhoSolid_;
    scalar rhoLiquid_;

    // Calculated Fields
    volScalarField Cp_;
    volScalarField rho_;
    volScalarField k_;

    void readDict();

public:

    //- Runtime type information
    TypeName("enthalpyPorosity");

    // Constructors
    pcmEnthalpyPorosityModel
    (
        const fvMesh& mesh,
        const solidThermo& thermo
    );

    // Destructor
    virtual ~pcmEnthalpyPorosityModel() = default;


    // Member Functions

    virtual bool active() const override
    {
        return active_;
    }

    virtual void correct() override;

    virtual tmp<volScalarField> CpEff() const override
    {
        return Cp_;
    }

    virtual tmp<volScalarField> rho() const override
    {
        return rho_;
    }

    virtual tmp<volScalarField> k() const override
    {
        return k_;
    }

    virtual tmp<volScalarField> latentHeatSource() const override;

    virtual tmp<volVectorField> momentumSource() const override;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
