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

#ifndef ehcPhaseChangeModel_H
#define ehcPhaseChangeModel_H

#include "phaseChangeModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class ehcPhaseChangeModel Declaration

 Description
    Effective Heat Capacity (EHC) model with quasi-Newton fixed-point secant latent heat
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
        const basicThermo& thermo
    );

    // Destructor
    virtual ~ehcPhaseChangeModel() = default;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
