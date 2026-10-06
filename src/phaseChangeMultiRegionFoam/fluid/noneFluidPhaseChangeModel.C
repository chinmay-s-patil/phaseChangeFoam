/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
\*---------------------------------------------------------------------------*/

#include "noneFluidPhaseChangeModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(noneFluidPhaseChangeModel, 0);
    addToRunTimeSelectionTable
    (
        fluidPhaseChangeModel,
        noneFluidPhaseChangeModel,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::noneFluidPhaseChangeModel::noneFluidPhaseChangeModel
(
    const fvMesh& mesh,
    const rhoReactionThermo& thermo,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    fluidPhaseChangeModel(mesh, thermo, U, phi)
{
    active_ = false;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::noneFluidPhaseChangeModel::~noneFluidPhaseChangeModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::noneFluidPhaseChangeModel::mDot() const
{
    return tmp<volScalarField>::New
    (
        IOobject
        (
            "mDot",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimMass/dimVolume/dimTime, Zero)
    );
}


Foam::tmp<Foam::volScalarField>
Foam::noneFluidPhaseChangeModel::divUSource() const
{
    return tmp<volScalarField>::New
    (
        IOobject
        (
            "divUSource",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimMass/dimVolume/dimTime, Zero)
    );
}


Foam::tmp<Foam::volScalarField>
Foam::noneFluidPhaseChangeModel::energySource() const
{
    return tmp<volScalarField>::New
    (
        IOobject
        (
            "energySource",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimEnergy/dimVolume/dimTime, Zero)
    );
}


void Foam::noneFluidPhaseChangeModel::correct()
{}


// ************************************************************************* //
