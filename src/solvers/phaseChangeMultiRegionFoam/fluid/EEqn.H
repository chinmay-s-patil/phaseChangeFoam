{
    volScalarField& he = thermo.he();

    if (pcmModel.active() && he.name() != "h")
    {
        FatalErrorInFunction
            << "PCM phase change model requires sensible enthalpy 'h' in thermo."
            << exit(FatalError);
    }

    tmp<volScalarField> tCpThermo = thermo.Cp();
    const volScalarField& CpThermo = tCpThermo();

    tmp<volScalarField> talphaEff =
        pcmModel.active()
      ? tmp<volScalarField>(pcmModel.k()/CpThermo + turbulence.alphat())
      : turbulence.alphaEff();
    const volScalarField& alphaEff = talphaEff();

    tmp<volScalarField> trhoEffCpRatio =
        pcmModel.active()
      ? tmp<volScalarField>(pcmModel.rhoEff()*pcmModel.CpEff()/CpThermo)
      : tmp<volScalarField>(rho);

    fvScalarMatrix EEqn
    (
        fvm::ddt(rho, he) + fvm::div(phi, he)
      + fvc::ddt(rho, K) + fvc::div(phi, K)
      + (
            he.name() == "e"
          ? fvc::div
            (
                fvc::absolute(phi/fvc::interpolate(rho), U),
                p,
                "div(phiv,p)"
            )
          : -dpdt
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
        EEqn.solve(he.select(finalIter));
        fvOptions.correct(he);

        thermo.correct();
        if (pcmModel.active())
        {
            pcmModel.postEnergySolve(he);
            thermo.correct();
        }
        rad.correct();

        Info<< "Min/max T:" << min(thermo.T()).value() << ' '
            << max(thermo.T()).value() << endl;
    }
}

