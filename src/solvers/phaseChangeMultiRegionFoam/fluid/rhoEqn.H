/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  www.openfoam.com
    \\  /    A nd           | Copyright (C) OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
\*---------------------------------------------------------------------------*/

{
    fvScalarMatrix rhoEqn
    (
        fvm::ddt(rho)
      + fvc::div(phi)
     ==
        fvOptions(rho)
    );

    if (phaseChange.active())
    {
        rhoEqn -= phaseChange.massSource();
    }

    fvOptions.constrain(rhoEqn);
    rhoEqn.solve();
    fvOptions.correct(rho);
}

// ************************************************************************* //
