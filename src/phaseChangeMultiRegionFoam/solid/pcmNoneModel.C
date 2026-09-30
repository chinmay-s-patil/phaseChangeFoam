/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Website:  www.openfoam.com
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

#include "pcmNoneModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(pcmNoneModel, 0);
    addToRunTimeSelectionTable(pcmPhaseChangeModel, pcmNoneModel, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pcmNoneModel::pcmNoneModel
(
    const fvMesh& mesh,
    const solidThermo& thermo
)
:
    pcmPhaseChangeModel(mesh, thermo)
{
    active_ = false;
}

// ************************************************************************* //
