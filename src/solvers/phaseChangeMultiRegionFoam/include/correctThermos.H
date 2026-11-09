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
    Modified 2026-10 by Chinmay S Patil: Updated fluid and solid thermophysical properties across all regions.

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

forAll(fluidRegions, i)
{
    rhoThermo& thermo = thermoFluid[i];
    radiation::radiationModel& rad = radiation[i];
    fv::options& fvOptions = fluidFvOptions[i];
    volScalarField& he = thermo.he();
    fvOptions.correct(he);
    thermo.correct();
    phaseChangeModel& pcmModel = pcmModelsFluid[i];
    if (pcmModel.active())
    {
        pcmModel.correct();
        pcmModel.postEnergySolve(he, finalIter);
        thermo.correct();
        pcmModel.applyKappaToThermo();
    }
    rad.correct();
}

forAll(solidRegions, i)
{
    solidThermo& thermo = thermos[i];
    fv::options& fvOptions = solidHeatSources[i];
    volScalarField& h = thermo.he();
    fvOptions.correct(h);
    thermo.correct();
    phaseChangeModel& pcModel = phaseChangeModels[i];
    if (pcModel.active())
    {
        pcModel.correct();
        pcModel.postEnergySolve(h, finalIter);
        thermo.correct();
        pcModel.applyKappaToThermo();
    }
}



