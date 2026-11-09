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
    Modified 2026-10 by Chinmay S Patil: Integrated steady-state phase change models and multi-region energy coupling.
    and multi-region energy coupling.

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

Application
    phaseChangeMultiRegionSimpleFoam

Description
    Steady-state multi-region solver for conjugate heat transfer with stateful
    generic phase-change (EHC / enthalpy-porosity) models.
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "turbulentFluidThermoModel.H"
#include "rhoReactionThermo.H"
#include "CombustionModel.H"
#include "fixedGradientFvPatchFields.H"
#include "regionProperties.H"
#include "compressibleCourantNo.H"
#include "solidRegionDiffNo.H"
#include "solidThermo.H"
#include "radiationModel.H"
#include "fvOptions.H"
#include "coordinateSystem.H"
#include "loopControl.H"
#include "pressureControl.H"
#include "phaseChangeModel.H"
#include "fluidPhaseChangeModel.H"
#include "leeFluidPhaseChangeModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Steady-state multi-region solver with phase-change heat transfer."
    );

    #define NO_CONTROL
    #define CREATE_MESH createMeshesPostProcess.H
    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMeshes.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"

    #include "createCoupledRegions.H"

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        const int oCorr = 0;
        const int nOuterCorr = 1;
        const bool finalIter = true;

        forAll(fluidRegions, i)
        {
            fvMesh& mesh = fluidRegions[i];

            Info<< "\nSolving for fluid region "
                << fluidRegions[i].name() << endl;
            #include "readFluidMultiRegionSIMPLEControls.H"
            #include "setRegionFluidFields.H"
            #include "solveFluid.H"
        }

        forAll(solidRegions, i)
        {
            fvMesh& mesh = solidRegions[i];

            #include "readSolidMultiRegionSIMPLEControls.H"
            #include "setRegionSolidFields.H"
            #include "solveSolid.H"
        }

        if (coupled)
        {
            Info<< "\nSolving energy coupled regions" << endl;
            fvMatrixAssemblyPtr->solve();
            #include "correctThermos.H"

            forAll(fluidRegions, i)
            {
                fvMesh& mesh = fluidRegions[i];

                #include "readFluidMultiRegionSIMPLEControls.H"
                #include "setRegionFluidFields.H"
                if (!frozenFlow)
                {
                    for (int corr=0; corr<nCorr; corr++)
                    {
                        #include "pEqn.H"
                    }
                    turbulence.correct();
                }
                rho = thermo.rho();
            }

            fvMatrixAssemblyPtr->clear();
        }

        // Additional loops for energy solution only
        {
            loopControl looping(runTime, "SIMPLE", "energyCoupling");

            while (looping.loop())
            {
                Info<< nl << looping << nl;

                forAll(fluidRegions, i)
                {
                    fvMesh& mesh = fluidRegions[i];

                    Info<< "\nSolving for fluid region "
                        << fluidRegions[i].name() << endl;
                    #include "readFluidMultiRegionSIMPLEControls.H"
                    #include "setRegionFluidFields.H"
                    frozenFlow = true;
                    #include "solveFluid.H"
                }

                forAll(solidRegions, i)
                {
                    fvMesh& mesh = solidRegions[i];

                    Info<< "\nSolving for solid region "
                        << solidRegions[i].name() << endl;
                    #include "readSolidMultiRegionSIMPLEControls.H"
                    #include "setRegionSolidFields.H"
                    #include "solveSolid.H"
                }

                if (coupled)
                {
                    Info<< "\nSolving energy coupled regions.. " << endl;
                    fvMatrixAssemblyPtr->solve();
                    #include "correctThermos.H"

                    forAll(fluidRegions, i)
                    {
                        #include "setRegionFluidFields.H"
                        rho = thermo.rho();
                    }

                    fvMatrixAssemblyPtr->clear();
                }
            }
        }

        // Re-synchronize phase-change models
        forAll(fluidRegions, i)
        {
            fluidPhaseChangeModel& phaseChange = fluidPhaseChangeModels[i];
            phaseChangeModel& pcmModel = pcmModelsFluid[i];
            if (phaseChange.active())
            {
                // Note: phaseChange.correct() for fluidPhaseChangeModel (e.g. Lee model) re-evaluates mDot using
                // final converged thermo state prior to runTime.write(), ensuring mDot is written consistently to disk.
                phaseChange.correct();
            }
            if (pcmModel.active())
            {
                pcmModel.correct();
                pcmModel.updateHistory();
            }
        }

        forAll(solidRegions, i)
        {
            phaseChangeModel& pcModel = phaseChangeModels[i];
            if (pcModel.active())
            {
                pcModel.correct();
                pcModel.updateHistory();
            }
        }


        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
