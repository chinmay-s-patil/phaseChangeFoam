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
    Modified 2026-10 by Chinmay S Patil: Populated fluid region meshes for steady solver.

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

    const wordList fluidNames(rp["fluid"]);

    PtrList<fvMesh> fluidRegions(fluidNames.size());

    forAll(fluidNames, i)
    {
        if (args.dryRun() || args.found("dry-run-write"))
        {
            Info
                << "Operating in 'dry-run' mode:"
                << " case will run for 1 time step.  "
                << "All checks assumed OK on a clean exit" << endl;

            FieldBase::allowConstructFromLargerSize = true;

            // Create a simplified 1D mesh and attempt to re-create boundary
            // conditions
            fluidRegions.set
            (
                i,
                new simplifiedMeshes::columnFvMesh(runTime, fluidNames[i])
            );

            // Stop after 1 iteration of the simplified mesh

            if (args.found("dry-run-write"))
            {
                // Using saWriteNow triggers function objects execute(), write()
                runTime.stopAt(Time::saWriteNow);

                // Make sure mesh gets output to the current time (since
                // instance no longer constant)
                fluidRegions[i].setInstance(runTime.timeName());
            }
            else
            {
                // Using saNoWriteNow triggers function objects execute(),
                // but not write()
                runTime.stopAt(Time::saNoWriteNow);
            }

            functionObject::outputPrefix = "postProcessing-dry-run";
        }
        else
        {
            Info<< "Create fluid mesh for region " << fluidNames[i]
                << " for time = " << runTime.timeName() << nl << endl;

            fluidRegions.set
            (
                i,
                new fvMesh
                (
                    IOobject
                    (
                        fluidNames[i],
                        runTime.timeName(),
                        runTime,
                        IOobject::MUST_READ
                    )
                )
            );
        }
    }
