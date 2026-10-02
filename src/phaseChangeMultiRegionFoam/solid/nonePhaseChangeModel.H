/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Website:  www.openfoam.com
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

#ifndef nonePhaseChangeModel_H
#define nonePhaseChangeModel_H

#include "phaseChangeModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class nonePhaseChangeModel Declaration

 Description
    Dummy phase-change model for solid regions with no phase change active.
    Registers baseline thermophysical fields (kEff, rhoEff, CpEff) to maintain
    compatibility with multi-region boundary condition lookups across regions.
\*---------------------------------------------------------------------------*/

class nonePhaseChangeModel
:
    public phaseChangeModel
{
public:

    //- Runtime type information
    TypeName("none");

private:
    volScalarField kEff_;
    volScalarField rhoEff_;
    volScalarField CpEff_;

public:

    // Constructors
    nonePhaseChangeModel
    (
        const fvMesh& mesh,
        const solidThermo& thermo
    );

    // Destructor
    virtual ~nonePhaseChangeModel() = default;


    // Member Functions

    virtual bool active() const override
    {
        return false;
    }

    virtual void correct() override
    {}

    virtual tmp<volScalarField> CpEff() const override
    {
        return CpEff_;
    }

    virtual tmp<volScalarField> rho() const override
    {
        return rhoEff_;
    }

    virtual tmp<volScalarField> k() const override
    {
        return kEff_;
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
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
