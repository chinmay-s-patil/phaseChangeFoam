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
