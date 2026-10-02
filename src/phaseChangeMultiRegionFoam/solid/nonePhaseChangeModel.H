/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Website:  www.openfoam.com
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

#ifndef pcmNoneModel_H
#define pcmNoneModel_H

#include "pcmPhaseChangeModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class pcmNoneModel Declaration
\*---------------------------------------------------------------------------*/

class pcmNoneModel
:
    public pcmPhaseChangeModel
{
public:

    //- Runtime type information
    TypeName("none");

    // Constructors
    pcmNoneModel
    (
        const fvMesh& mesh,
        const solidThermo& thermo
    );

    // Destructor
    virtual ~pcmNoneModel() = default;


    // Member Functions

    virtual bool active() const override
    {
        return false;
    }

    virtual void correct() override
    {}

    virtual tmp<volScalarField> CpEff() const override
    {
        return thermo_.Cp();
    }

    virtual tmp<volScalarField> rho() const override
    {
        return thermo_.rho();
    }

    virtual tmp<volScalarField> k() const override
    {
        return thermo_.kappa();
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
