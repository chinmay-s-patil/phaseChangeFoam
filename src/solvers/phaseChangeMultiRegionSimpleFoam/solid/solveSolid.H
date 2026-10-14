{
    if (pcModel.active())
    {
        pcModel.correct();
    }

    tmp<volScalarField> tkEff = pcModel.active() ? pcModel.k() : tmp<volScalarField>(thermo.kappa());
    tmp<volScalarField> tCpThermo = thermo.Cp();
    tmp<volScalarField> trhoEff = pcModel.active() ? pcModel.rhoEff() : tmp<volScalarField>(thermo.rho());
    volScalarField alphaEff("alphaEff", tkEff() / tCpThermo());

    for (int nonOrth=0; nonOrth<=nNonOrthCorr; ++nonOrth)
    {
        fvScalarMatrix hEqn
        (
            (
                thermo.isotropic()
              ? -fvm::laplacian(betav*alphaEff, h, "laplacian(alpha,h)")
              : -fvm::laplacian(betav*taniAlpha(), h, "laplacian(alpha,h)")
            )
          ==
            fvOptions(trhoEff(), h)
        );

        if (pcModel.active())
        {
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
            Info<< "\nSolving for solid region "
                << solidRegions[i].name() << endl;

            hEqn.solve();
            fvOptions.correct(h);
            thermo.correct();

            if (pcModel.active())
            {
                pcModel.correct();
                pcModel.postEnergySolve(h);
                thermo.correct();
                pcModel.applyKappaToThermo();
            }

            Info<< "Min/max T:" << min(thermo.T()).value() << ' '
                << max(thermo.T()).value() << endl;
        }
    }
}

