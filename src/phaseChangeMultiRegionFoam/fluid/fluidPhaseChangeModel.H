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
        fluidPhaseChangeModel(const fluidPhaseChangeModel&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const fluidPhaseChangeModel&) = delete;


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

            //- Net mixture mass source rate per unit volume [kg/(m^3 s)]
            //  Placed in mass continuity / pressure equation (pEqn.H: - phaseChange.massSource()).
            //  Zero for internal single-phase species conversion (conserves total domain mass);
            //  non-zero if liquid enters from an external source (film, pool, or spray).
            virtual tmp<volScalarField> massSource() const = 0;

            //- Species mass source rate per unit volume by specie index [kg/(m^3 s)]
            //  Placed on RHS of species transport equation (YEqn.H: + phaseChange.speciesSource(i)).
            //  Positive for species creation (e.g. vapor created by evaporation),
            //  negative for species consumption (e.g. liquid species depleted in cell).
            virtual tmp<volScalarField> speciesSource(const label specieIndex) const = 0;

            //- Species mass source rate per unit volume by specie name [kg/(m^3 s)]
            //  Default implementation maps specieName to index in thermo composition.
            virtual tmp<volScalarField> speciesSource(const word& specieName) const;

            //- Latent and mass enthalpy source/sink term for energy equation [W/m^3]
            //  Placed on RHS of energy equation (EEqn.H: EEqn == ... + phaseChange.energySource()).
            //  Note: Evaporation absorbs latent heat, so latent term is negative (- mDot_evap * L).
            //  If net mixture mass enters the region (massSource() > 0), energySource() also
            //  includes the enthalpy carried by added mass (+ mDot_added * h_added).
            virtual tmp<volScalarField> energySource() const = 0;

            //- Update phase change state and fields
            virtual void correct() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
