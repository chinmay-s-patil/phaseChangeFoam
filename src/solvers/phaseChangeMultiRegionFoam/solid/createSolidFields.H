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
    Modified 2026-10 by Chinmay S Patil: Instantiated solid-region phaseChangeModel pointer lists and alpha fields.

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

    // Initialise solid field pointer lists
    PtrList<coordinateSystem> coordinates(solidRegions.size());
    PtrList<solidThermo> thermos(solidRegions.size());
    PtrList<radiation::radiationModel> radiations(solidRegions.size());
    PtrList<fv::options> solidHeatSources(solidRegions.size());
    PtrList<volScalarField> betavSolid(solidRegions.size());
    PtrList<volSymmTensorField> aniAlphas(solidRegions.size());
    PtrList<phaseChangeModel> phaseChangeModels(solidRegions.size());


    // Populate solid field pointer lists
    forAll(solidRegions, i)
    {
        Info<< "*** Reading solid mesh thermophysical properties for region "
            << solidRegions[i].name() << nl << endl;

        Info<< "    Adding to thermos\n" << endl;
        thermos.set(i, solidThermo::New(solidRegions[i]));

        Info<< "    Adding phase change model\n" << endl;
        phaseChangeModels.set(i, phaseChangeModel::New(solidRegions[i], thermos[i]));

        if (phaseChangeModels[i].active() && !thermos[i].isotropic())
        {
            FatalErrorInFunction
                << "Phase change is active in solid region '" << solidRegions[i].name()
                << "', but the region thermophysical model is anisotropic (!thermo.isotropic()).\n"
                << "Phase change models compute isotropic effective thermal conductivity k(alpha);\n"
                << "anisotropic phase change heat conduction is not supported."
                << exit(FatalError);
        }

        Info<< "    Adding to radiations\n" << endl;
        radiations.set(i, radiation::radiationModel::New(thermos[i].T()));

        Info<< "    Adding fvOptions\n" << endl;
        solidHeatSources.set
        (
            i,
            new fv::options(solidRegions[i])
        );

        if (!thermos[i].isotropic())
        {
            Info<< "    Adding coordinateSystems\n" << endl;
            coordinates.set
            (
                i,
                coordinateSystem::New
                (
                    solidRegions[i],
                    thermos[i],
                    coordinateSystem::typeName
                )
            );

            tmp<volVectorField> tkappaByCp =
                thermos[i].Kappa()/thermos[i].Cp();

            aniAlphas.set
            (
                i,
                new volSymmTensorField
                (
                    IOobject
                    (
                        "Anialpha",
                        runTime.timeName(),
                        solidRegions[i],
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    solidRegions[i],
                    dimensionedSymmTensor(tkappaByCp().dimensions(), Zero),
                    fvPatchFieldBase::zeroGradientType()
                )
            );

            aniAlphas[i].primitiveFieldRef() =
                coordinates[i].transformPrincipal
                (
                    solidRegions[i].cellCentres(),
                    tkappaByCp()
                );
            aniAlphas[i].correctBoundaryConditions();

        }

        IOobject betavSolidIO
        (
            "betavSolid",
            runTime.timeName(),
            solidRegions[i],
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        );

        if (betavSolidIO.typeHeaderOk<volScalarField>(true))
        {
            betavSolid.set
            (
                i,
                new volScalarField(betavSolidIO, solidRegions[i])
            );
        }
        else
        {
            betavSolid.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        "betavSolid",
                        runTime.timeName(),
                        solidRegions[i],
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    solidRegions[i],
                    dimensionedScalar("1", dimless, scalar(1))
                )
            );
        }
    }
