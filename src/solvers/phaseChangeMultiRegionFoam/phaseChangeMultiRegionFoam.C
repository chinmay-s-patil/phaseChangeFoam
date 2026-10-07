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
    phaseChangeMultiRegionFoam

Description
    Multi-region solver for conjugate heat transfer with stateful
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
        "Multi-region solver with phase-change heat transfer."
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
    #include "createTimeControls.H"
    #include "readSolidTimeControls.H"
    #include "compressibleMultiRegionCourantNo.H"
    #include "solidRegionDiffusionNo.H"
    #include "setInitialMultiRegionDeltaT.H"

    #include "createCoupledRegions.H"

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "readSolidTimeControls.H"
        #include "readPIMPLEControls.H"

        #include "compressibleMultiRegionCourantNo.H"
        #include "solidRegionDiffusionNo.H"
        #include "setMultiRegionDeltaT.H"

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        if (nOuterCorr != 1)
        {
            forAll(fluidRegions, i)
            {
                #include "storeOldFluidFields.H"
            }
        }

        PtrList<volScalarField> T_outer_prev_solid(solidRegions.size());
        PtrList<volScalarField> T_outer_prev_fluid(fluidRegions.size());
        if (coupled)
        {
            forAll(solidRegions, i)
            {
                T_outer_prev_solid.set
                (
                    i,
                    new volScalarField
                    (
                        IOobject("T_outer_prev_solid", runTime.timeName(), solidRegions[i]),
                        thermos[i].T()
                    )
                );
            }
            forAll(fluidRegions, i)
            {
                T_outer_prev_fluid.set
                (
                    i,
                    new volScalarField
                    (
                        IOobject("T_outer_prev_fluid", runTime.timeName(), fluidRegions[i]),
                        thermoFluid[i].T()
                    )
                );
            }
        }

        // --- PIMPLE loop
        for (int oCorr=0; oCorr<nOuterCorr; ++oCorr)
        {
            const bool finalIter = (oCorr == nOuterCorr-1);

            forAll(fluidRegions, i)
            {
                fvMesh& mesh = fluidRegions[i];

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
                Info<< "\nSolving energy coupled regions " << endl;
                fvMatrixAssemblyPtr->solve();
                #include "correctThermos.H"

                scalar outerTol = pimple.getOrDefault<scalar>("outerTolerance", 1e-4);

                scalar maxDeltaT_outer = 0.0;
                forAll(solidRegions, i)
                {
                    maxDeltaT_outer = max(maxDeltaT_outer, max(mag(thermos[i].T() - T_outer_prev_solid[i])).value());
                    T_outer_prev_solid[i] = thermos[i].T();
                }
                forAll(fluidRegions, i)
                {
                    maxDeltaT_outer = max(maxDeltaT_outer, max(mag(thermoFluid[i].T() - T_outer_prev_fluid[i])).value());
                    T_outer_prev_fluid[i] = thermoFluid[i].T();
                }

                Info<< "Coupled outer convergence max DeltaT = " << maxDeltaT_outer << " K (tol = " << outerTol << " K)" << endl;

                if (oCorr > 0 && maxDeltaT_outer < outerTol)
                {
                    Info<< "Coupled outer iterations CONVERGED at outer corrector " << (oCorr + 1) << endl;
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
                    break;
                }

                if (finalIter && nOuterCorr > 1 && maxDeltaT_outer >= outerTol)
                {
                    WarningInFunction
                        << "Coupled outer iterations did not converge within " << nOuterCorr
                        << " outer correctors (residual max DeltaT = " << maxDeltaT_outer
                        << " K > tol = " << outerTol << " K)" << endl;
                }

                forAll(fluidRegions, i)
                {
                    fvMesh& mesh = fluidRegions[i];

                    #include "readFluidMultiRegionPIMPLEControls.H"
                    #include "setRegionFluidFields.H"
                    if (!frozenFlow)
                    {
                        Info<< "\nSolving for fluid region "
                            << fluidRegions[i].name() << endl;
                        // --- PISO loop
                        for (int corr=0; corr<nCorr; corr++)
                        {
                            #include "pEqn.H"
                        }
                        turbulence.correct();
                    }

                    rho = thermo.rho();
                    Info<< "Min/max T:" << min(thermo.T()).value() << ' '
                        << max(thermo.T()).value() << endl;
                }

                fvMatrixAssemblyPtr->clear();
            }

            // Additional loops for energy solution only
            if (coupled && !oCorr && nOuterCorr > 1)
            {
                loopControl looping(runTime, pimple, "energyCoupling");

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
                        Info<< "\nSolving energy coupled regions " << endl;
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
        }

        // Re-synchronize phase-change models with final converged T^n before committing history.
        forAll(fluidRegions, i)
        {
            fluidPhaseChangeModel& phaseChange = fluidPhaseChangeModels[i];
            phaseChangeModel& pcmModel = pcmModelsFluid[i];
            if (phaseChange.active())
            {
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
