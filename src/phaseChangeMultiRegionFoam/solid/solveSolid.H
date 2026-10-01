{
    int nNonLinearCorr = pcmModel.active() ? pimple.getOrDefault<int>("nNonLinearCorrectors", 3) : 1;

    for (int nlCorr = 0; nlCorr < nNonLinearCorr; ++nlCorr)
    {
        if (pcmModel.active())
        {
            pcmModel.correct();
        }

        tmp<volScalarField> tCpThermo = thermo.Cp();
        tmp<volScalarField> tCpEff = pcmModel.active() ? pcmModel.CpEff() : tmp<volScalarField>(thermo.Cp());
        tmp<volScalarField> trhoPCM = pcmModel.active() ? pcmModel.rho() : tmp<volScalarField>(thermo.rho());
        tmp<volScalarField> tkPCM = pcmModel.active() ? pcmModel.k() : tmp<volScalarField>(thermo.kappa());

        // Thermal diffusivity alpha = k / Cp_sensible in kg/(m s) for enthalpy-based equation
        volScalarField alpha("alpha", tkPCM() / tCpThermo());

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
            break;
        }
        else
        {
            if (nlCorr == 0)
            {
                Info<< "\nSolving for solid region "<< solidRegions[i].name() << endl;
            }

            if (finalIter && (nlCorr == nNonLinearCorr - 1))
            {
                mesh.data().setFinalIteration(true);
            }

            hEqn.solve(h.select(finalIter));

            fvOptions.correct(h);

            thermo.correct();

            if (finalIter && (nlCorr == nNonLinearCorr - 1))
            {
                mesh.data().setFinalIteration(false);
            }
        }
    }

    if (!coupled)
    {
        Info<< "Min/max T:" << min(thermo.T()).value() << ' '
            << max(thermo.T()).value() << endl;
    }
}
