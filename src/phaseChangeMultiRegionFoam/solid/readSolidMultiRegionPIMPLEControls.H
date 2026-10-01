const dictionary& pimple = mesh.solutionDict().subDict("PIMPLE");

int nNonOrthCorr =
    pimple.getOrDefault<int>("nNonOrthogonalCorrectors", 0);

int nNonLinearCorr =
    pimple.getOrDefault<int>("nNonLinearCorrectors", 3);
