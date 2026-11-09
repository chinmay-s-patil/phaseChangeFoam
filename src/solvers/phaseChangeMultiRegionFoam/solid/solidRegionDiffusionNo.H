/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2019,2022 OpenCFD Ltd.
    Copyright (C) 2026 Chinmay S Patil
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM, with modifications for phaseChangeFoam.
    SPDX-License-Identifier: GPL-3.0-or-later
    Modified 2026-10 by Chinmay S Patil: Evaluated solid region diffusion number using phase-dependent thermal diffusivity.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

scalar DiNum = -GREAT;

forAll(solidRegions, i)
{
    //- Note: do not use setRegionSolidFields.H to avoid double registering Cp
    //#include "setRegionSolidFields.H"
    const solidThermo& thermo = thermos[i];
    const phaseChangeModel& pcModel = phaseChangeModels[i];

    tmp<volScalarField> magKappa =
        pcModel.active()
      ? pcModel.k()
      : (thermo.isotropic() ? thermo.kappa() : tmp<volScalarField>(mag(thermo.Kappa())));

    tmp<volScalarField> tcp = pcModel.active() ? pcModel.CpEff() : tmp<volScalarField>(thermo.Cp());
    tmp<volScalarField> trho = pcModel.active() ? pcModel.rhoEff() : tmp<volScalarField>(thermo.rho());

    volScalarField Cprho("Cprho", trho() * tcp());

    DiNum = max
    (
        solidRegionDiffNo
        (
            solidRegions[i],
            runTime,
            Cprho,
            magKappa()
        ),
        DiNum
    );

}
