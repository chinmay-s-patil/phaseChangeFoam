{
    if (pcModel.active())
    {
        pcModel.correct();
    }

    tmp<volScalarField> tkEff = pcModel.active() ? pcModel.k() : tmp<volScalarField>(thermo.kappa());
    tmp<volScalarField> tCpThermo = thermo.Cp();
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
            fvOptions(rho, h)
        );

        if (pcModel.active())
        {
            hEqn += pcModel.latentHeatSource()
                  + fvm::Sp(pcModel.latentHeatSp(), h)
                  - pcModel.latentHeatSp() * h;
        }

        hEqn.relax();

        fvOptions.constrain(hEqn);

        if (coupled)
        {
            fvMatrixAssemblyPtr->addFvMatrix(hEqn);
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
            }

            Info<< "Min/max T:" << min(thermo.T()).value() << ' '
                << max(thermo.T()).value() << endl;
        }
    }
}

