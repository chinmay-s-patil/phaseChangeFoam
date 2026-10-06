/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Class
    Foam::fluidPhaseChangeModel

Description
    Abstract base class for fluid-region phase-change and evaporation models.
    Provides mass source rate mDot(), velocity divergence source divUSource(),
    and energy-equation latent heat source/sink terms.

\*---------------------------------------------------------------------------*/

#ifndef fluidPhaseChangeModel_H
#define fluidPhaseChangeModel_H

#include "IOdictionary.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "rhoReactionThermo.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class fluidPhaseChangeModel Declaration
\*---------------------------------------------------------------------------*/

class fluidPhaseChangeModel
:
    public IOdictionary
{
protected:

    // Protected data

        const fvMesh& mesh_;
        const rhoReactionThermo& thermo_;
        const volVectorField& U_;
        const surfaceScalarField& phi_;

        dictionary dict_;
        bool active_;


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        fluidPhaseChangeModel(const fluidPhaseChangeModel&);

        //- Disallow default bitwise assignment
        void operator=(const fluidPhaseChangeModel&);


public:

    //- Runtime type information
    TypeName("fluidPhaseChangeModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            fluidPhaseChangeModel,
            dictionary,
            (
                const fvMesh& mesh,
                const rhoReactionThermo& thermo,
                const volVectorField& U,
                const surfaceScalarField& phi
            ),
            (mesh, thermo, U, phi)
        );


    // Constructors

        //- Construct from mesh, thermo, U, phi
        fluidPhaseChangeModel
        (
            const fvMesh& mesh,
            const rhoReactionThermo& thermo,
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    // Selector

        static autoPtr<fluidPhaseChangeModel> New
        (
            const fvMesh& mesh,
            const rhoReactionThermo& thermo,
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~fluidPhaseChangeModel();


    // Member Functions

        // Access

            //- Return mesh reference
            const fvMesh& mesh() const
            {
                return mesh_;
            }

            //- Return true if model is active
            virtual bool active() const
            {
                return active_;
            }


        // Abstract interface methods

            //- Mass transfer rate per unit volume [kg/(m^3 s)]
            virtual tmp<volScalarField> mDot() const = 0;

            //- Mass divergence / volume expansion source [kg/(m^3 s)]
            virtual tmp<volScalarField> divUSource() const = 0;

            //- Latent heat source/sink term for energy equation [W/m^3]
            virtual tmp<volScalarField> energySource() const = 0;

            //- Update phase change state and fields
            virtual void correct() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
