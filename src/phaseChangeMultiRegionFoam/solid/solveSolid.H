{
    if (pcmModel.active())
    {
        pcmModel.correct();
    }

    tmp<volScalarField> tCpThermo = thermo.Cp();
    tmp<volScalarField> tCpEff = pcmModel.active() ? pcmModel.CpEff() : tmp<volScalarField>(thermo.Cp());
    tmp<volScalarField> trhoPCM = pcmModel.active() ? pcmModel.rho() : tmp<volScalarField>(thermo.rho());
    tmp<volScalarField> tkPCM = pcmModel.active() ? pcmModel.k() : tmp<volScalarField>(thermo.kappa());

    // Sensible thermal diffusivity alpha = k / (rho * Cp_sensible)
    volScalarField alpha("alpha", tkPCM() / (trhoPCM() * tCpThermo()));

    // Effective heat capacity scaling ratio for transient term: rho * (CpEff / Cp_sensible)
    volScalarField rhoCpRatio("rhoCpRatio", trhoPCM() * (tCpEff() / tCpThermo()));

    fvScalarMatrix hEqn
    (
        fvm::ddt(betav*rhoCpRatio, h)
        - (
            thermo.isotropic()
            ? fvm::laplacian(betav*alpha, h, "laplacian(alpha,h)")
            : fvm::laplacian(betav*taniAlpha(), h, "laplacian(alpha,h)")
        )
        ==
        fvOptions(trhoPCM(), h)
    );

    if (pcmModel.active())
    {
        hEqn += pcmModel.latentHeatSource();
    }

    hEqn.relax();

    fvOptions.constrain(hEqn);

    if (coupled)
    {
        fvMatrixAssemblyPtr->addFvMatrix(hEqn);
    }
    else
    {
        Info<< "\nSolving for solid region "<< solidRegions[i].name() << endl;

        if (finalIter)
        {
            mesh.data().setFinalIteration(true);
        }

        hEqn.solve(h.select(finalIter));

        fvOptions.correct(h);

        thermo.correct();

        Info<< "Min/max T:" << min(thermo.T()).value() << ' '
            << max(thermo.T()).value() << endl;

        if (finalIter)
        {
            mesh.data().setFinalIteration(false);
        }
    }
}
