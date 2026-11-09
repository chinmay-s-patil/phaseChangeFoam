/*--------------------------------*- C++ -*----------------------------------*\
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
        const basicThermo& thermo
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
