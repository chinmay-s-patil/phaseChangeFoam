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

#include "ehcPhaseChangeModel.H"
#include "calculatedFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(ehcPhaseChangeModel, 0);
    addToRunTimeSelectionTable(phaseChangeModel, ehcPhaseChangeModel, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ehcPhaseChangeModel::ehcPhaseChangeModel
(
    const fvMesh& mesh,
    const basicThermo& thermo
)
:
    phaseChangeModel(mesh, thermo, true)
{
    if (mesh.foundObject<volVectorField>("U") && !suppressConvection_)
    {
        FatalErrorInFunction
            << "Model 'EHC' (ehcPhaseChangeModel) does not include momentum Darcy drag."
            << "\nConvection cannot be enabled on a fluid region with EHC."
            << "\nUse model 'EHCDarcy' or 'enthalpyPorosity', or set 'convection { suppress true; }'."
            << exit(FatalError);
    }


    // Validate that all non-constraint patches are 'calculated'
    auto checkCalculatedPatches = [&](const volScalarField& f)
    {
        forAll(f.boundaryField(), pI)
        {
            const fvPatchScalarField& pf = f.boundaryField()[pI];
            if
            (
                pf.type() != calculatedFvPatchScalarField::typeName
             && !polyPatch::constraintType(pf.patch().patch().type())
            )
            {
                FatalErrorInFunction
                    << "Field '" << f.name() << "' patch '" << pf.patch().name()
                    << "' has type '" << pf.type() << "' but must be 'calculated'.\n"
                    << "Remove or correct the patch entry in the restart file."
                    << exit(FatalError);
            }
        }
    };
    checkCalculatedPatches(phaseFraction_);
    checkCalculatedPatches(heatingTrajectory_);
    checkCalculatedPatches(T_reversal_);
    checkCalculatedPatches(phaseState_);
    checkCalculatedPatches(Cp_);
    checkCalculatedPatches(rho_);
    checkCalculatedPatches(k_);
}

// ************************************************************************* //
