/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Class
    Foam::constantSourceFluidPhaseChangeModel

Description
    Fluid phase-change test model returning constant, user-specified volumetric
    mass, species, and energy sources for plumbing verification.

\*---------------------------------------------------------------------------*/

#ifndef constantSourceFluidPhaseChangeModel_H
#define constantSourceFluidPhaseChangeModel_H

#include "fluidPhaseChangeModel.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class constantSourceFluidPhaseChangeModel Declaration
\*---------------------------------------------------------------------------*/

class constantSourceFluidPhaseChangeModel
:
    public fluidPhaseChangeModel
{

private:

    // Private Data

        word speciesName_;
        label specieIndex_;

        volScalarField massSource_;
        volScalarField energySource_;
        volScalarField speciesSource_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        constantSourceFluidPhaseChangeModel(const constantSourceFluidPhaseChangeModel&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const constantSourceFluidPhaseChangeModel&) = delete;


public:

    //- Runtime type information
    TypeName("constantSource");


    // Constructors

        //- Construct from mesh, thermo, U, phi
        constantSourceFluidPhaseChangeModel
        (
            const fvMesh& mesh,
            const rhoReactionThermo& thermo,
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~constantSourceFluidPhaseChangeModel();


    // Member Functions

        // Access

            //- Always return true for active model
            virtual bool active() const override
            {
                return active_;
            }


        // Abstract interface implementations

            //- Return net mixture mass source rate [kg/(m^3 s)]
            virtual tmp<volScalarField> massSource() const override;

            //- Return species mass source rate by specie index [kg/(m^3 s)]
            virtual tmp<volScalarField> speciesSource(const label specieIndex) const override;

            //- Return energy source [W/m^3]
            virtual tmp<volScalarField> energySource() const override;

            //- Update phase change state and fields
            virtual void correct() override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
