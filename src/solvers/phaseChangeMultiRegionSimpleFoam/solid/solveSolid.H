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
    Modified 2026-10 by Chinmay S Patil: Implemented steady solid energy equation with pcmModel latent heat source terms.

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
                pcModel.postEnergySolve(h);
                thermo.correct();
                pcModel.applyKappaToThermo();
            }

            Info<< "Min/max T:" << min(thermo.T()).value() << ' '
                << max(thermo.T()).value() << endl;
        }
    }
}

