/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 Chinmay S Patil
-------------------------------------------------------------------------------
License
    This file is part of phaseChangeFoam.
    SPDX-License-Identifier: GPL-3.0-or-later

    phaseChangeFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    phaseChangeFoam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
    or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with phaseChangeFoam.  If not, see <http://www.gnu.org/licenses/>.
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
            dict_.lookupOrDefault<scalar>
            (
                "massSource",
                dict_.lookupOrDefault<scalar>("speciesSource", 0.0)
            )
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
            dict_.lookupOrDefault<scalar>
            (
                "speciesSource",
                dict_.lookupOrDefault<scalar>("massSource", 0.0)
            )
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
