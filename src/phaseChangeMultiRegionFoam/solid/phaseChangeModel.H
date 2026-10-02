/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Website:  www.openfoam.com
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

#ifndef phaseChangeModel_H
#define phaseChangeModel_H

#include "fvCFD.H"
#include "solidThermo.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class phaseChangeModel Declaration

 Description
    Base class for general multi-region phase-change thermal models
    (e.g. solidification/melting, alloys, phase-change materials,
    ice-water transition, and general thermal energy storage/transfer).

    Provides standardized interfaces for effective thermophysical properties,
    latent heat source terms, and phase state machine transitions.
\*---------------------------------------------------------------------------*/

class phaseChangeModel
{
protected:

    // Protected data
    const fvMesh& mesh_;
    const solidThermo& thermo_;

    bool active_;

    // Fields
    volScalarField phaseFraction_; // Liquid/transformed phase fraction alpha in [0, 1]
    volScalarField phaseState_;    // 0 = solid, 1 = melting/transitioning, 2 = liquid, 3 = freezing

    //- True if phaseFraction was found on disk (or restored from legacy liquidFraction)
    bool phaseFractionRestored_;

    //- Convection suppression flag (dict: phaseChange.convection.suppress, default true).
    //  Regions solved by this application carry no momentum equation, so convection
    //  is always suppressed; false is rejected at read time.
    bool suppressConvection_;

    //- Read optional 'convection { suppress <bool>; }' sub-dictionary
    void readConvectionDict(const dictionary& pcDict);

public:

    //- Runtime type information
    TypeName("phaseChangeModel");

    //- Declare runtime construction table
    declareRunTimeSelectionTable
    (
        autoPtr,
        phaseChangeModel,
        dictionary,
        (
            const fvMesh& mesh,
            const solidThermo& thermo
        ),
        (mesh, thermo)
    );

    // Constructors
    phaseChangeModel
    (
        const fvMesh& mesh,
        const solidThermo& thermo
    );

    // Selector
    static autoPtr<phaseChangeModel> New
    (
        const fvMesh& mesh,
        const solidThermo& thermo
    );

    // Destructor
    virtual ~phaseChangeModel() = default;


    // Member Functions

    //- Check if phase change is active for this region
    virtual bool active() const
    {
        return active_;
    }

    //- Correct phase fraction, state machine, and thermo fields
    virtual void correct() = 0;

    //- Update trajectory history at completed timestep boundaries
    virtual void updateHistory() {}

    //- Whether this model suppresses fluid convection (solid-only phase change)
    virtual bool suppressConvection() const
    {
        return suppressConvection_;
    }

    //- Access phase fraction field
    virtual const volScalarField& phaseFraction() const
    {
        return phaseFraction_;
    }

    //- Backward compatibility alias for liquid fraction
    virtual const volScalarField& liquidFraction() const
    {
        return phaseFraction_;
    }

    //- Access phase state field
    virtual const volScalarField& phaseState() const
    {
        return phaseState_;
    }

    //- Access effective heat capacity field
    virtual tmp<volScalarField> CpEff() const = 0;

    //- Access density field
    virtual tmp<volScalarField> rho() const = 0;

    //- Access thermal conductivity field
    virtual tmp<volScalarField> k() const = 0;

    //- Access latent heat source term
    virtual tmp<volScalarField> latentHeatSource() const = 0;

    //- Access momentum source term (default: zero, override in models with Darcy damping)
    virtual tmp<volVectorField> momentumSource() const
    {
        return tmp<volVectorField>::New
        (
            IOobject("momentumSource", mesh_.time().timeName(), mesh_),
            mesh_,
            dimensionedVector("zero", dimForce/dimVolume, Zero)
        );
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
