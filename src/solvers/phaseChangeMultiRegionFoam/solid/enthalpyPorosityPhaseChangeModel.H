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
    Voller–Prakash Enthalpy-Porosity phase change model.
    Uses source-based latent energy advection L*[ddt(rho,f) + div(phi,f)],
    implicit Darcy momentum drag fvm::Sp, and a post-energy-solve Newton
    liquid fraction update with enthalpy correction h -= L*df.
\*---------------------------------------------------------------------------*/

class enthalpyPorosityPhaseChangeModel
:
    public phaseChangeModel
{
private:

    // Parameters
    scalar lambda_; // Newton relaxation factor
    scalar L_;      // Latent heat of fusion [J/kg]
    scalar Ts_;     // Solidus temperature [K]
    scalar Tl_;     // Liquidus temperature [K]
    scalar Cu_;     // Darcy penalty coefficient [kg/(m3 s)]
    scalar q_;      // Smoothing constant
    scalar Cp_;     // Single specific heat capacity [J/(kg K)]
    scalar kSolid_; // Solid thermal conductivity [W/(m K)]
    scalar kLiquid_;// Liquid thermal conductivity [W/(m K)]

    void readEPDict();

public:

    //- Runtime type information
    TypeName("enthalpyPorosity");

    // Constructors
    enthalpyPorosityPhaseChangeModel
    (
        const fvMesh& mesh,
        const basicThermo& thermo
    );

    // Destructor
    virtual ~enthalpyPorosityPhaseChangeModel() = default;

    // Member Functions

    //- Correct phase change model fields
    virtual void correct() override;

    //- Update history (no-op for EP, OpenFOAM handles oldTime level)
    virtual void updateHistory() override {}

    //- Latent heat source term: L * [ fvc::ddt(rho, f) + fvc::div(phi, f) ]
    virtual tmp<volScalarField> latentHeatSource() const override;

    //- Implicit latent heat linearisation coefficient (zero for EP)
    virtual tmp<volScalarField> latentHeatSp() const override;

    //- Access implicit Darcy momentum damping coefficient [kg/(m^3 s)]
    virtual tmp<volScalarField> momentumSp() const override;

    //- Access density field for energy equation (returns background thermo.rho() for EP)
    virtual tmp<volScalarField> rhoEff() const override;

    //- Access effective heat capacity field (returns background thermo.Cp() for EP)
    virtual tmp<volScalarField> CpEff() const override;

    //- Access thermal conductivity field k(f) = f*kl + (1-f)*ks
    virtual tmp<volScalarField> k() const override;

    //- Post-energy-solve Newton liquid fraction update and enthalpy adjustment
    virtual void postEnergySolve(volScalarField& he) override;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
