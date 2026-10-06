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
    const dictionary& pimple = solutionDict.subDict("PIMPLE");
    const int nOuterCorr = pimple.getOrDefault<int>("nOuterCorrectors", 1);
    
    bool anyPhaseChangeActive = false;
    forAll(solidRegions, i)
    {
        if (phaseChangeModels[i].active())
        {
            anyPhaseChangeActive = true;
            if (coupled)
            {
                Info<< "Phase change active in coupled solid region "
                    << solidRegions[i].name()
                    << " with nOuterCorrectors = " << nOuterCorr << endl;
            }
        }
    }
    
    const bool isMultiRegion = (fluidRegions.size() + solidRegions.size()) > 1;
    if (anyPhaseChangeActive && (coupled || isMultiRegion) && nOuterCorr < 2)
    {
        FatalErrorInFunction
            << "Phase change is active in a multi-region coupled domain, but nOuterCorrectors = "
            << nOuterCorr << " < 2.\n"
            << "Multi-region phase change requires nOuterCorrectors >= 2 "
            << "to converge non-linear latent heat and temperature across region boundaries."
            << exit(FatalError);
    }
}
