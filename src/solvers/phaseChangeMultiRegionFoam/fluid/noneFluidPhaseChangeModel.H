/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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
