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
    Modified 2026-10 by Chinmay S Patil: Integrated phase change mass source (rhoPhCh) into pressure correction equation.

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

bool closedVolume = p_rgh.needReference();
dimensionedScalar compressibility = fvc::domainIntegrate(psi);
bool compressible = (compressibility.value() > SMALL);

rho = thermo.rho();

// Thermodynamic density needs to be updated by psi*d(p) after the
// pressure solution
const volScalarField psip0(psi*p);

volScalarField rAU("rAU", 1.0/UEqn.A());
surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho*rAU, "interpolate(rhorAUf)"));
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p_rgh));

surfaceScalarField phig(-rhorAUf*ghf*fvc::snGrad(rho)*mesh.magSf());

surfaceScalarField phiHbyA
(
    "phiHbyA",
    (
        fvc::flux(rho*HbyA)
      + MRF.zeroFilter(rhorAUf*fvc::ddtCorr(rho, U, phi))
    )
  + phig
);

MRF.makeRelative(fvc::interpolate(rho), phiHbyA);

// Update the pressure BCs to ensure flux consistency
constrainPressure(p_rgh, rho, U, phiHbyA, rhorAUf, MRF);

{
    fvScalarMatrix p_rghDDtEqn
    (
        fvc::ddt(rho) + psi*correction(fvm::ddt(p_rgh))
      + fvc::div(phiHbyA)
    );

    if (phaseChange.active())
    {
        p_rghDDtEqn -= phaseChange.massSource();
    }

    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix p_rghEqn
        (
            p_rghDDtEqn
          - fvm::laplacian(rhorAUf, p_rgh)
        );

        p_rghEqn.setReference
        (
            pRefCell,
            compressible ? getRefCellValue(p_rgh, pRefCell) : pRefValue
        );

        p_rghEqn.solve
        (
            p_rgh.select
            (
                (
                    oCorr == nOuterCorr-1
                 && corr == nCorr-1
                 && nonOrth == nNonOrthCorr
                )
            )
        );

        if (nonOrth == nNonOrthCorr)
        {
            phi = phiHbyA + p_rghEqn.flux();

            p_rgh.relax();

            U = HbyA
              + rAU*fvc::reconstruct((phig + p_rghEqn.flux())/rhorAUf);
            U.correctBoundaryConditions();
            fvOptions.correct(U);
            K = 0.5*magSqr(U);
        }
    }

    p = p_rgh + rho*gh;

}

pressureControl.limit(p);

// For closed-volume cases adjust the pressure and density levels
// to obey overall mass continuity
if (closedVolume)
{
    if (!compressible)
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
    }
    else
    {
        p += (initialMass - fvc::domainIntegrate(psi*p))
            /compressibility;
        thermo.correctRho(psi*p - psip0, rhoMin, rhoMax);
        rho = thermo.rho();
        p_rgh = p - rho*gh;
        p_rgh.correctBoundaryConditions();
    }
}
else
{
    thermo.correctRho(psi*p - psip0,  rhoMin, rhoMax);
}

#include "rhoEqn.H"
#include "compressibleContinuityErrors.H"

rho = thermo.rho();

// Update pressure time derivative if needed
if (thermo.dpdt())
{
    dpdt = fvc::ddt(p);
}
