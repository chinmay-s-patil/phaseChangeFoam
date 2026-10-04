/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Website:  www.openfoam.com
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

#ifndef enthalpyPorosityPhaseChangeModel_H
#define enthalpyPorosityPhaseChangeModel_H

#include "phaseChangeModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class enthalpyPorosityPhaseChangeModel Declaration

 Description
    Enthalpy-Porosity model for general phase-change simulations with exact
    Newton/secant linearization, stateful hysteresis, path-integral thermophysics,
    and Darcy momentum damping force for fluid regions.
\*---------------------------------------------------------------------------*/

class enthalpyPorosityPhaseChangeModel
:
    public phaseChangeModel
{
private:

    // Parameters - Darcy & Buoyancy
    scalar beta_; // Thermal expansion coefficient [1/K]
    scalar Cu_;   // Darcy penalty coefficient [kg/(m3 s)]
    scalar q_;    // Smoothing constant

    void readPorosityDict();

public:

    //- Runtime type information
    TypeName("enthalpyPorosity");

    // Constructors
    enthalpyPorosityPhaseChangeModel
    (
        const fvMesh& mesh,
        const solidThermo& thermo
    );

    // Destructor
    virtual ~enthalpyPorosityPhaseChangeModel() = default;


    // Member Functions

    virtual tmp<volVectorField> momentumSource() const override;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
