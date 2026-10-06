/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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
            #include "readFluidMultiRegionPIMPLEControls.H"
            #include "setRegionFluidFields.H"
            #include "solveFluid.H"
        }

        forAll(solidRegions, i)
        {
            fvMesh& mesh = solidRegions[i];

            #include "readSolidMultiRegionPIMPLEControls.H"
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

                #include "readFluidMultiRegionPIMPLEControls.H"
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
                    #include "readFluidMultiRegionPIMPLEControls.H"
                    #include "setRegionFluidFields.H"
                    frozenFlow = true;
                    #include "solveFluid.H"
                }

                forAll(solidRegions, i)
                {
                    fvMesh& mesh = solidRegions[i];

                    Info<< "\nSolving for solid region "
                        << solidRegions[i].name() << endl;
                    #include "readSolidMultiRegionPIMPLEControls.H"
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
