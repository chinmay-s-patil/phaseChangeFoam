/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 Chinmay S Patil
-------------------------------------------------------------------------------
License
    This file is part of phaseChangeFoam.
    SPDX-License-Identifier: GPL-3.0-or-later

    phaseChangeFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    phaseChangeFoam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
    or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with phaseChangeFoam.  If not, see <http://www.gnu.org/licenses/>.
Class
    Foam::noneFluidPhaseChangeModel

Description
    Dummy / inactive fluid phase-change model returning zero sources.

\*---------------------------------------------------------------------------*/

#ifndef noneFluidPhaseChangeModel_H
#define noneFluidPhaseChangeModel_H

#include "fluidPhaseChangeModel.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class noneFluidPhaseChangeModel Declaration
\*---------------------------------------------------------------------------*/

class noneFluidPhaseChangeModel
:
    public fluidPhaseChangeModel
{

private:

    // Private Data

        volScalarField zeroMassSource_;
        volScalarField zeroEnergySource_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        noneFluidPhaseChangeModel(const noneFluidPhaseChangeModel&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const noneFluidPhaseChangeModel&) = delete;


public:

    //- Runtime type information
    TypeName("none");


    // Constructors

        //- Construct from mesh, thermo, U, phi
        noneFluidPhaseChangeModel
        (
            const fvMesh& mesh,
            const rhoReactionThermo& thermo,
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~noneFluidPhaseChangeModel();


    // Member Functions

        // Access

            //- Always return false for inactive model
            virtual bool active() const override
            {
                return false;
            }


        // Abstract interface implementations

            //- Return zero mixture mass source rate [kg/(m^3 s)]
            virtual tmp<volScalarField> massSource() const override;

            //- Return zero species mass source rate by index [kg/(m^3 s)]
            virtual tmp<volScalarField> speciesSource(const label specieIndex) const override;

            //- Return zero energy source [W/m^3]
            virtual tmp<volScalarField> energySource() const override;

            //- No-op update
            virtual void correct() override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
