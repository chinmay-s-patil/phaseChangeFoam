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
