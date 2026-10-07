{
    volScalarField& he = thermo.he();

    tmp<volScalarField> tCpThermo = thermo.Cp();
    const volScalarField& CpThermo = tCpThermo();

    tmp<volScalarField> talphaEff =
        pcmModel.active()
      ? tmp<volScalarField>(volScalarField("alphaEffPCM", pcmModel.k()/CpThermo + turb.alphat()))
      : turb.alphaEff();
    const volScalarField& alphaEff = talphaEff();

    fvScalarMatrix EEqn
    (
        fvm::div(phi, he)
      + (
            he.name() == "e"
          ? fvc::div(phi, volScalarField("Ekp", 0.5*magSqr(U) + p/rho))
          : fvc::div(phi, volScalarField("K", 0.5*magSqr(U)))
        )
      - fvm::laplacian(alphaEff, he)
     ==
        rho*(U&g)
      + rad.Sh(thermo, he)
      + Qdot
      + fvOptions(rho, he)
    );



    if (phaseChange.active())
    {
        EEqn -= phaseChange.energySource();
    }

    if (pcmModel.active())
    {
        EEqn += pcmModel.latentHeatSource()
              + fvm::Sp(pcmModel.latentHeatSp(), he)
              - pcmModel.latentHeatSp() * he;
    }

    EEqn.relax();


    fvOptions.constrain(EEqn);

    if (coupled)
    {
        fvMatrixAssemblyPtr->addFvMatrix(EEqn);
    }
    else
    {
        EEqn.solve();

        fvOptions.correct(he);

        thermo.correct();
        rad.correct();

        Info<< "Min/max T:" << min(thermo.T()).value() << ' '
            << max(thermo.T()).value() << endl;
    }
}
