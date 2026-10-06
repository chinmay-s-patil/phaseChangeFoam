/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Class
    Foam::leeFluidPhaseChangeModel

Description
    Lee volumetric phase-change model for fluid regions.
    Uses mixture density scaling (rho * Y_specie) for single-phase in-region
    species conversion and supports optional pressure-dependent Tsat(p).

    Physics Caveats & Model Scope:
    1. Single-Region Species Conversion Scope (massSource_ = 0):
       massSource() returns zero because this model converts species within a single phase/region
       (conserving total cell mass). Inter-region mass exchange (e.g. mass influx from a liquid
       film or pool region into the gas region) is not implemented in this Lee model and
       requires an inter-region phase change model.
    2. Boiling-Type Model Assumption:
       Tsat(p) is evaluated against total cell pressure p. Phase change occurs
       when mixture T > Tsat(p). Evaporation into non-condensable gas below the
       boiling point (e.g. water evaporating into dry air) is driven by vapor
       partial pressure (p_v,sat(T) vs X_v * p), which requires a boundary or
       inter-region mass transfer model rather than total-pressure Lee.
    3. Latent Heat vs Specie Enthalpy of Formation (Hf):
       The energy source term energySource() = - mDot * L assumes sensible
       enthalpy thermo where latent heat L is explicitly provided. If species
       thermodynamic data already encodes phase enthalpy difference in Hf
       (e.g., Hf_v - Hf_l = L), adding - mDot * L can double-count latent heat.

\*---------------------------------------------------------------------------*/

#ifndef leeFluidPhaseChangeModel_H
#define leeFluidPhaseChangeModel_H

#include "fluidPhaseChangeModel.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class leeFluidPhaseChangeModel Declaration
\*---------------------------------------------------------------------------*/

class leeFluidPhaseChangeModel
:
    public fluidPhaseChangeModel
{

private:

    // Private Data

        word liquidName_;
        word vaporName_;
        label liquidIndex_;
        label vaporIndex_;

        dimensionedScalar C_evap_;
        dimensionedScalar C_cond_;
        dimensionedScalar latentHeat_;
        dimensionedScalar TsatRef_;
        dimensionedScalar pRef_;
        bool enableTsatP_;

        volScalarField mDot_;
        volScalarField massSource_;
        volScalarField energySource_;
        volScalarField liquidSource_;
        volScalarField vaporSource_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        leeFluidPhaseChangeModel(const leeFluidPhaseChangeModel&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const leeFluidPhaseChangeModel&) = delete;


public:

    //- Runtime type information
    TypeName("Lee");


    // Constructors

        //- Construct from mesh, thermo, U, phi
        leeFluidPhaseChangeModel
        (
            const fvMesh& mesh,
            const rhoReactionThermo& thermo,
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~leeFluidPhaseChangeModel();


    // Member Functions

        // Access

            //- Return true if active
            virtual bool active() const override
            {
                return active_;
            }

            //- Return phase change mass rate field [kg/(m^3 s)]
            const volScalarField& mDot() const
            {
                return mDot_;
            }


        // Abstract interface implementations

            //- Return net mixture mass source rate [kg/(m^3 s)].
            //  Returns zero for in-region species conversion (net cell mass conserved).
            //  Inter-region mass influx is not implemented in this model.
            virtual tmp<volScalarField> massSource() const override;

            //- Return species mass source rate by specie index [kg/(m^3 s)]
            virtual tmp<volScalarField> speciesSource(const label specieIndex) const override;

            //- Return energy source [W/m^3]
            virtual tmp<volScalarField> energySource() const override;

            //- Update phase change state and source fields
            virtual void correct() override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
