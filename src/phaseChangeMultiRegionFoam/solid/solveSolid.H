{
    int nNonLinearCorr = pimple.getOrDefault<int>("nNonLinearCorrectors", 25);
    scalar nonLinearTol = pimple.getOrDefault<scalar>("nonLinearTolerance", 1e-4);

    volScalarField Tprev("Tprev", thermo.T());

    int actualIter = 0;
    scalar maxDeltaT = 0.0;
    for (int nlCorr = 0; nlCorr < nNonLinearCorr; ++nlCorr)
    {
        actualIter = nlCorr + 1;

        if (pcModel.active())
        {
            pcModel.correct();
        }

        tmp<volScalarField> tCpThermo = thermo.Cp();
        tmp<volScalarField> tCpEff = pcModel.active() ? pcModel.CpEff() : tmp<volScalarField>(thermo.Cp());
        tmp<volScalarField> trhoEff = pcModel.active() ? pcModel.rhoEff() : tmp<volScalarField>(thermo.rho());
        tmp<volScalarField> tkEff = pcModel.active() ? pcModel.k() : tmp<volScalarField>(thermo.kappa());

        // Thermal diffusivity alpha = k / Cp_thermo so that alpha * grad(h) == k * grad(T) exactly
        volScalarField alpha("alpha", tkEff() / tCpThermo());

        // Effective thermal capacity storage ratio (rho_eff * Cp_eff) / Cp_thermo
        volScalarField rhoCpRatio("rhoCpRatio", trhoEff() * tCpEff() / tCpThermo());

        fvScalarMatrix hEqn
        (
            fvm::ddt(betav*rhoCpRatio, h)
            - (
                thermo.isotropic()
                ? fvm::laplacian(betav*alpha, h, "laplacian(alpha,h)")
                : fvm::laplacian(betav*taniAlpha(), h, "laplacian(alpha,h)")
            )
            ==
            fvOptions(trhoEff(), h)
        );

        if (pcModel.active())
        {
            // Linearized latent heat source: Su_phys + fvm::Sp(Sp, h) - Sp * h_k
            // At convergence (h = h_k), fvm::Sp(Sp, h) and -Sp * h_k cancel out,
            // leaving exact Su_phys = -rho * L * dAlpha/dt on RHS.
            hEqn += pcModel.latentHeatSource()
                  + fvm::Sp(pcModel.latentHeatSp(), h)
                  - pcModel.latentHeatSp() * h;
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

            if (pcModel.active())
            {
                pcModel.correct();
            }
            else
            {
                thermo.correct();
            }

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
