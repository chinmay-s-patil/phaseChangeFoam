/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Website:  www.openfoam.com
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

#ifndef ehcDarcyPhaseChangeModel_H
#define ehcDarcyPhaseChangeModel_H

#include "phaseChangeModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class ehcDarcyPhaseChangeModel Declaration

 Description
    EHCDarcy model combining path-integral Effective Heat Capacity (EHC) for
    energy phase change with implicit Darcy momentum damping force for fluid regions.
\*---------------------------------------------------------------------------*/

class ehcDarcyPhaseChangeModel
:
    public phaseChangeModel
{
private:

    // Parameters - Darcy
    scalar Cu_;   // Darcy penalty coefficient [kg/(m3 s)]
    scalar q_;    // Smoothing constant


    void readPorosityDict();

public:

    //- Runtime type information
    TypeName("EHCDarcy");

    // Constructors
    ehcDarcyPhaseChangeModel
    (
        const fvMesh& mesh,
        const basicThermo& thermo
    );

    // Destructor
    virtual ~ehcDarcyPhaseChangeModel() = default;

    // Member Functions

    //- Access implicit Darcy momentum damping coefficient [kg/(m^3 s)]
    virtual tmp<volScalarField> momentumSp() const override;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
