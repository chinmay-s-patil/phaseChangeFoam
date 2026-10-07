/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Website:  www.openfoam.com
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

#ifndef enthalpyPorosityPhaseChangeModel_H
#define enthalpyPorosityPhaseChangeModel_H

#include "ehcDarcyPhaseChangeModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class enthalpyPorosityPhaseChangeModel Declaration

 Description
    Enthalpy-Porosity model alias deriving from EHCDarcy for backward
    compatibility in case dictionaries.
\*---------------------------------------------------------------------------*/

class enthalpyPorosityPhaseChangeModel
:
    public ehcDarcyPhaseChangeModel
{
public:

    //- Runtime type information
    TypeName("enthalpyPorosity");

    // Constructors
    enthalpyPorosityPhaseChangeModel
    (
        const fvMesh& mesh,
        const basicThermo& thermo
    )
    :
        ehcDarcyPhaseChangeModel(mesh, thermo)
    {}

    // Destructor
    virtual ~enthalpyPorosityPhaseChangeModel() = default;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
