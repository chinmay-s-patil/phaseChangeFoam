/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2019,2022 OpenCFD Ltd.
    Copyright (C) 2026 Chinmay S Patil
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM, with modifications for phaseChangeFoam.
    SPDX-License-Identifier: GPL-3.0-or-later
    Modified 2026-10 by Chinmay S Patil: Added phase change sensible enthalpy source terms (mDotL, Sh) and linearisation.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

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

    bool compressible = (fvc::domainIntegrate(thermo.psi()).value() > SMALL);

    fvScalarMatrix EEqn
    (
        fvm::ddt(rho, he) + fvm::div(phi, he)
      - fvm::laplacian(alphaEff, he)
     ==
        rad.Sh(thermo, he)
      + Qdot
      + fvOptions(rho, he)
    );

    if (compressible)
    {
        EEqn +=
            fvc::ddt(rho, K) + fvc::div(phi, K)
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
          - rho*(U&g);
    }

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
            pcmModel.postEnergySolve(he, finalIter);
            thermo.correct();
            pcmModel.applyKappaToThermo();
        }
        rad.correct();

        Info<< "Min/max T:" << min(thermo.T()).value() << ' '
            << max(thermo.T()).value() << endl;
    }
}

