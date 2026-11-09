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
