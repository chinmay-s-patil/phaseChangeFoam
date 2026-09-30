/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Website:  www.openfoam.com
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

#ifndef pcmPhaseChangeModel_H
#define pcmPhaseChangeModel_H

#include "fvCFD.H"
#include "solidThermo.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class pcmPhaseChangeModel Declaration
\*---------------------------------------------------------------------------*/

class pcmPhaseChangeModel
{
protected:

    // Protected data
    const fvMesh& mesh_;
    const solidThermo& thermo_;

    bool active_;

    // Fields
    volScalarField liquidFraction_; // Phase fraction alpha_L in [0, 1]
    volScalarField phaseState_;     // 0 = solid, 1 = melting, 2 = liquid, 3 = freezing

public:

    //- Runtime type information
    TypeName("pcmPhaseChangeModel");

    //- Declare runtime construction table
    declareRunTimeSelectionTable
    (
        autoPtr,
        pcmPhaseChangeModel,
        dictionary,
        (
            const fvMesh& mesh,
            const solidThermo& thermo
        ),
        (mesh, thermo)
    );

    // Constructors
    pcmPhaseChangeModel
    (
        const fvMesh& mesh,
        const solidThermo& thermo
    );

    // Selector
    static autoPtr<pcmPhaseChangeModel> New
    (
        const fvMesh& mesh,
        const solidThermo& thermo
    );

    // Destructor
    virtual ~pcmPhaseChangeModel() = default;


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

    //- Access liquid fraction field
    virtual const volScalarField& liquidFraction() const
    {
        return liquidFraction_;
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

    //- Access momentum source term
    virtual tmp<volVectorField> momentumSource() const = 0;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
