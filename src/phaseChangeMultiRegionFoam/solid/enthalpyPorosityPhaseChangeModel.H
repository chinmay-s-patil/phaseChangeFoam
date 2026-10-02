/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Website:  www.openfoam.com
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

#ifndef enthalpyPorosityPhaseChangeModel_H
#define enthalpyPorosityPhaseChangeModel_H

#include "phaseChangeModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class enthalpyPorosityPhaseChangeModel Declaration

 Description
    Enthalpy-Porosity model for general phase-change simulations (solidification,
    melting, alloy casting, and phase-change thermal storage).
    Note: Full classical enthalpy-porosity with momentum damping (Darcy source)
    requires a fluid region where UEqn is solved. In solid phase-change regions,
    sensible heat conduction and latent heat source terms are active.
\*---------------------------------------------------------------------------*/

class enthalpyPorosityPhaseChangeModel
:
    public phaseChangeModel
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
    scalar relax_; // Non-linear relaxation for phase fraction update

    // Parameters - Density & Thermo
    word densityModel_;
    scalar rhoRef_;
    scalar rhoSolid_;
    scalar rhoLiquid_;

    // Calculated Fields
    volScalarField Cp_;
    volScalarField rho_;
    volScalarField k_;

    //- Phase fraction at start of step (explicit; avoids oldTime() capture ordering)
    volScalarField phaseFraction_old_;

    void readDict();

public:

    //- Runtime type information
    TypeName("enthalpyPorosity");

    // Constructors
    enthalpyPorosityPhaseChangeModel
    (
        const fvMesh& mesh,
        const solidThermo& thermo
    );

    // Destructor
    virtual ~enthalpyPorosityPhaseChangeModel() = default;


    // Member Functions

    virtual bool active() const override
    {
        return active_;
    }

    virtual void correct() override;

    virtual void updateHistory() override;

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
