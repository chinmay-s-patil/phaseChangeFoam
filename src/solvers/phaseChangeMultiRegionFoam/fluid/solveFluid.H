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
    Modified 2026-10 by Chinmay S Patil: Invoked phaseChange.correct() and pcmModel.correct() within the fluid solution step.

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

if (finalIter)
{
    mesh.data().setFinalIteration(true);
}

if (oCorr == 0)
{
    bool compressible = (fvc::domainIntegrate(thermo.psi()).value() > SMALL);
    if (!frozenFlow && compressible)
    {
        #include "rhoEqn.H"
    }
    phaseChange.correct();
    phaseChange.checkMassConservation();
}

if (pcmModel.active())
{
    pcmModel.correct();
}

if (frozenFlow)
{
    #include "EEqn.H"
}
else
{
    #include "UEqn.H"
    #include "YEqn.H"
    #include "EEqn.H"

    if (!coupled)
    {
        Info<< "\nSolving for fluid region " << fluidRegions[i].name() << endl;

        // --- PISO loop
        for (int corr=0; corr<nCorr; corr++)
        {
            #include "pEqn.H"
        }

        turbulence.correct();

        rho = thermo.rho();
    }
}

if (finalIter)
{
    mesh.data().setFinalIteration(false);
}
