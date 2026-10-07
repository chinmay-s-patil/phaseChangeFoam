{
    int nNonLinearCorr = pimple.getOrDefault<int>("nNonLinearCorrectors", 25);
    scalar nonLinearTol = pimple.getOrDefault<scalar>("nonLinearTolerance", 1e-4);
    scalar nonLinearTolAlpha = pimple.getOrDefault<scalar>("nonLinearToleranceAlpha", 1e-4);

    volScalarField Tprev("Tprev", thermo.T());
    volScalarField alphaPrev
    (
        IOobject
        (
            "alphaPrev",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0)
    );
    if (pcModel.active())
    {
        alphaPrev = pcModel.phaseFraction();
    }

    int actualIter = 0;
    scalar maxDeltaT = 0.0;
    scalar maxDeltaAlpha = 0.0;
    for (int nlCorr = 0; nlCorr < nNonLinearCorr; ++nlCorr)
    {
        actualIter = nlCorr + 1;

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
            betav*rhoCpRatio*fvm::ddt(h)
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
            // Scaled by solid volume fraction betav
            hEqn += betav * pcModel.latentHeatSource()
                  + fvm::Sp(betav * pcModel.latentHeatSp(), h)
                  - betav * pcModel.latentHeatSp() * h;
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
                pcModel.postEnergySolve(h);
                thermo.correct();
            }
            else
            {
                thermo.correct();
            }


            maxDeltaT = max(mag(thermo.T() - Tprev)).value();
            Tprev = thermo.T();

            if (pcModel.active())
            {
                maxDeltaAlpha = max(mag(pcModel.phaseFraction() - alphaPrev)).value();
                alphaPrev = pcModel.phaseFraction();
            }

            if (finalIter && (nlCorr == nNonLinearCorr - 1))
            {
                mesh.data().setFinalIteration(false);
            }

            bool converged = (maxDeltaT < nonLinearTol);
            if (pcModel.active())
            {
                converged = converged && (maxDeltaAlpha < nonLinearTolAlpha);
            }

            if (converged && nlCorr > 0)
            {
                break;
            }
        }
    }

    if (!coupled)
    {
        bool converged = (maxDeltaT < nonLinearTol);
        if (pcModel.active())
        {
            converged = converged && (maxDeltaAlpha < nonLinearTolAlpha);
        }

        if (actualIter == nNonLinearCorr && !converged)
        {
            if (pcModel.active())
            {
                WarningInFunction
                    << "Solid region " << solidRegions[i].name()
                    << " did not converge within " << nNonLinearCorr
                    << " non-linear iterations (maxDeltaT = " << maxDeltaT
                    << " K > tol = " << nonLinearTol << " K, maxDeltaAlpha = "
                    << maxDeltaAlpha << " > tol = " << nonLinearTolAlpha << ")" << endl;
            }
            else
            {
                WarningInFunction
                    << "Solid region " << solidRegions[i].name()
                    << " did not converge within " << nNonLinearCorr
                    << " non-linear iterations (residual maxDeltaT = " << maxDeltaT
                    << " K > tol = " << nonLinearTol << " K)" << endl;
            }
        }

        Info<< "Min/max T:" << min(thermo.T()).value() << ' '
            << max(thermo.T()).value()
            << " (non-linear iters: " << actualIter << ")" << endl;
    }
}
