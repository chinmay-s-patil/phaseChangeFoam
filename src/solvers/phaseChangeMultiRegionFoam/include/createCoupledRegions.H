bool coupled = false;
{
    fvSolution solutionDict(runTime);
    solutionDict.readIfPresent("coupledEnergyField", coupled);
}

autoPtr<fvMatrix<scalar>> fvMatrixAssemblyPtr;

forAll(fluidRegions, i)
{
    const rhoThermo& thermo = refCast<const rhoThermo>(thermoFluid[i]);
    const auto& bpsi = thermo.T().boundaryField();

    forAll(bpsi, patchI)
    {
        if (bpsi[patchI].useImplicit())
        {
            coupled = true;
        }
    }
}

forAll(solidRegions, i)
{
    solidThermo& thermo = thermos[i];
    const auto& bpsi = thermo.T().boundaryField();

    forAll(bpsi, patchI)
    {
        if (bpsi[patchI].useImplicit())
        {
            coupled = true;
        }
    }
}

forAll(fluidRegions, i)
{
    const rhoThermo& thermo = refCast<const rhoThermo>(thermoFluid[i]);
    if (coupled)
    {
        Info<< "Create fvMatrixAssembly." << endl;
        fvMatrixAssemblyPtr.reset
        (
            new fvMatrix<scalar>
            (
                thermo.he(),
                dimEnergy/dimTime
            )
        );
        break;
    }
}

if (coupled && !fvMatrixAssemblyPtr)
{
    forAll(solidRegions, i)
    {
        solidThermo& thermo = thermos[i];
        Info<< "Create fvMatrixAssembly." << endl;
        fvMatrixAssemblyPtr.reset
        (
            new fvMatrix<scalar>
            (
                thermo.he(),
                dimEnergy/dimTime
            )
        );
        break;
    }
}

{
    fvSolution solutionDict(runTime);
    const dictionary* pimplePtr =
        solutionDict.found("PIMPLE")
      ? &solutionDict.subDict("PIMPLE")
      : (solutionDict.found("SIMPLE") ? &solutionDict.subDict("SIMPLE") : nullptr);

    const int nOuterCorr =
        pimplePtr
      ? pimplePtr->getOrDefault<int>("nOuterCorrectors", 1)
      : 1;

    bool anySolidPhaseChangeActive = false;
    forAll(solidRegions, i)
    {
        if (phaseChangeModels[i].active())
        {
            anySolidPhaseChangeActive = true;
            if (coupled)
            {
                Info<< "Phase change active in coupled solid region "
                    << solidRegions[i].name()
                    << " with nOuterCorrectors = " << nOuterCorr << endl;
            }
        }
    }

    bool anyFluidPhaseChangeActive = false;
    forAll(fluidRegions, i)
    {
        if (pcmModelsFluid[i].active() || fluidPhaseChangeModels[i].active())
        {
            anyFluidPhaseChangeActive = true;
            Info<< "Phase change active in fluid region "
                << fluidRegions[i].name()
                << " with nOuterCorrectors = " << nOuterCorr << endl;
        }
    }

    const bool anyPhaseChangeActive = anySolidPhaseChangeActive || anyFluidPhaseChangeActive;
    const bool isMultiRegion = (fluidRegions.size() + solidRegions.size()) > 1;

    if (anyPhaseChangeActive && (anyFluidPhaseChangeActive || coupled || isMultiRegion) && nOuterCorr < 2)
    {
        FatalErrorInFunction
            << "Phase change is active in a fluid, coupled, or multi-region domain, but nOuterCorrectors = "
            << nOuterCorr << " < 2.\n"
            << "Fluid, coupled, and multi-region phase change require nOuterCorrectors >= 2 "
            << "to non-linearly iterate and converge latent heat release and temperature within each time step."
            << exit(FatalError);
    }
}
