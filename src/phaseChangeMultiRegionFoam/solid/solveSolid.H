{
    int nNonLinearCorr = pimple.getOrDefault<int>("nNonLinearCorrectors", 25);
    scalar nonLinearTol = pimple.getOrDefault<scalar>("nonLinearTolerance", 1e-4);

    volScalarField Tprev("Tprev", thermo.T());

    int actualIter = 0;
    scalar maxDeltaT = 0.0;
    for (int nlCorr = 0; nlCorr < nNonLinearCorr; ++nlCorr)
    {
        actualIter = nlCorr + 1;

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
            if (pcmModel.active() && (nlCorr == 0))
            {
                WarningInFunction
                    << "Solid region " << solidRegions[i].name()
                    << " has active phase change in a coupled multi-region domain. "
                    << "Single-pass matrix assembly relies on outer PIMPLE iterations (nOuterCorrectors > 1) "
                    << "to converge phase change non-linearity across regions." << endl;
            }
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

            maxDeltaT = max(mag(thermo.T() - Tprev)).value();
            Tprev = thermo.T();

            if (finalIter && (nlCorr == nNonLinearCorr - 1))
            {
                mesh.data().setFinalIteration(false);
            }

            if (maxDeltaT < nonLinearTol && nlCorr > 0)
            {
                break;
            }
        }
    }

    if (!coupled)
    {
        if (actualIter == nNonLinearCorr && maxDeltaT >= nonLinearTol)
        {
            WarningInFunction
                << "Solid region " << solidRegions[i].name()
                << " did not converge within " << nNonLinearCorr
                << " non-linear iterations (residual maxDeltaT = " << maxDeltaT
                << " K > tol = " << nonLinearTol << " K)" << endl;
        }

        Info<< "Min/max T:" << min(thermo.T()).value() << ' '
            << max(thermo.T()).value()
            << " (non-linear iters: " << actualIter << ")" << endl;
    }
}
