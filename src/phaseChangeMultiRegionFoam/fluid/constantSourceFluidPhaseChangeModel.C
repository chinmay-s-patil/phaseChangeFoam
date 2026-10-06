/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
\*---------------------------------------------------------------------------*/

#include "constantSourceFluidPhaseChangeModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(constantSourceFluidPhaseChangeModel, 0);
    addToRunTimeSelectionTable
    (
        fluidPhaseChangeModel,
        constantSourceFluidPhaseChangeModel,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::constantSourceFluidPhaseChangeModel::constantSourceFluidPhaseChangeModel
(
    const fvMesh& mesh,
    const rhoReactionThermo& thermo,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    fluidPhaseChangeModel(mesh, thermo, U, phi),
    speciesName_(dict_.getOrDefault<word>("speciesName", "H2O")),
    specieIndex_(thermo_.composition().species().find(speciesName_)),
    massSource_
    (
        IOobject
        (
            "constantMassSource",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "massSourceVal",
            dimMass/dimVolume/dimTime,
            dict_.getOrDefault<scalar>("massSource", 0.0)
        )
    ),
    energySource_
    (
        IOobject
        (
            "constantEnergySource",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "energySourceVal",
            dimEnergy/dimVolume/dimTime,
            dict_.getOrDefault<scalar>("energySource", -1e6)
        )
    ),
    speciesSource_
    (
        IOobject
        (
            "constantSpeciesSource",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "speciesSourceVal",
            dimMass/dimVolume/dimTime,
            dict_.getOrDefault<scalar>("speciesSource", 1e-3)
        )
    )
{
    active_ = true;
    if (specieIndex_ == -1)
    {
        Info<< "constantSourceFluidPhaseChangeModel: specie " << speciesName_
            << " not found in thermo composition. speciesSource will be zero." << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::constantSourceFluidPhaseChangeModel::~constantSourceFluidPhaseChangeModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::constantSourceFluidPhaseChangeModel::massSource() const
{
    return massSource_;
}


Foam::tmp<Foam::volScalarField>
Foam::constantSourceFluidPhaseChangeModel::speciesSource(const label specieIndex) const
{
    if (specieIndex == specieIndex_ && specieIndex_ != -1)
    {
        return speciesSource_;
    }
    else
    {
        return zeroSpeciesSource_;
    }
}


Foam::tmp<Foam::volScalarField>
Foam::constantSourceFluidPhaseChangeModel::energySource() const
{
    return energySource_;
}


void Foam::constantSourceFluidPhaseChangeModel::correct()
{}


// ************************************************************************* //
