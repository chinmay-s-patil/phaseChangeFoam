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

    fvOptions.constrain(rhoEqn);
    rhoEqn.solve();
    fvOptions.correct(rho);
}

// ************************************************************************* //
