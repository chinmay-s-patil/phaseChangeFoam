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
    fluidPhaseChangeModel(mesh, thermo, U, phi),
    zeroMassSource_
    (
        IOobject
        (
            "zeroMassSource",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimMass/dimVolume/dimTime, Zero)
    ),
    zeroEnergySource_
    (
        IOobject
        (
            "zeroEnergySource",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimEnergy/dimVolume/dimTime, Zero)
    )
{
    active_ = false;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::noneFluidPhaseChangeModel::~noneFluidPhaseChangeModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::noneFluidPhaseChangeModel::massSource() const
{
    return zeroMassSource_;
}


Foam::tmp<Foam::volScalarField>
Foam::noneFluidPhaseChangeModel::speciesSource(const label specieIndex) const
{
    return zeroMassSource_;
}


Foam::tmp<Foam::volScalarField>
Foam::noneFluidPhaseChangeModel::energySource() const
{
    return zeroEnergySource_;
}


void Foam::noneFluidPhaseChangeModel::correct()
{}


// ************************************************************************* //
