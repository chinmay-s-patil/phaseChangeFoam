/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Website:  www.openfoam.com
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

#ifndef ehcPhaseChangeModel_H
#define ehcPhaseChangeModel_H

#include "phaseChangeModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class ehcPhaseChangeModel Declaration

 Description
    Effective Heat Capacity (EHC) model with exact Newton/secant latent heat
    linearisation, stateful hysteresis, and path-integral thermophysics
    for solid phase-change heat transfer.
\*---------------------------------------------------------------------------*/

class ehcPhaseChangeModel
:
    public phaseChangeModel
{
public:

    //- Runtime type information
    TypeName("EHC");

    // Constructors
    ehcPhaseChangeModel
    (
        const fvMesh& mesh,
        const solidThermo& thermo
    );

    // Destructor
    virtual ~ehcPhaseChangeModel() = default;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
