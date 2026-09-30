/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Website:  www.openfoam.com
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

#include "pcmEnthalpyPorosityModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(pcmEnthalpyPorosityModel, 0);
    addToRunTimeSelectionTable(pcmPhaseChangeModel, pcmEnthalpyPorosityModel, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pcmEnthalpyPorosityModel::pcmEnthalpyPorosityModel
(
    const fvMesh& mesh,
    const solidThermo& thermo
)
:
    pcmPhaseChangeModel(mesh, thermo),
    Tlm_(303.15),
    Tum_(313.15),
    Lm_(163000.0),
    beta_(50e-6),
    Cu_(1e5),
    q_(1e-2),
    rhoRef_(1967.0),
    rhoSolid_(1967.0),
    rhoLiquid_(1850.0),
    Cp_
    (
        IOobject("CpEP", mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
        thermo.Cp()
    ),
    rho_
    (
        IOobject("rhoEP", mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
        thermo.rho()
    ),
    k_
    (
        IOobject("kEP", mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
        thermo.kappa()
    )
{
    active_ = true;
    readDict();
    correct();
}

void Foam::pcmEnthalpyPorosityModel::readDict()
{
    IOobject dictIO
    (
        "phaseChangeDict",
        mesh_.time().constant(),
        mesh_,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    IOdictionary phaseChangeDict(dictIO);
    const dictionary& pcDict = phaseChangeDict.subDict("phaseChange");

    if (pcDict.found("melting"))
    {
        const dictionary& meltDict = pcDict.subDict("melting");
        Tlm_ = meltDict.lookupOrDefault<scalar>("T_lowerBound", 303.15);
        Tum_ = meltDict.lookupOrDefault<scalar>("T_upperBound", 313.15);
        Lm_ = meltDict.lookupOrDefault<scalar>("latentHeat", 163000.0);
    }

    if (pcDict.found("enthalpyPorosity"))
    {
        const dictionary& epDict = pcDict.subDict("enthalpyPorosity");
        beta_ = epDict.lookupOrDefault<scalar>("beta", 50e-6);
        Cu_ = epDict.lookupOrDefault<scalar>("Cu", 1e5);
        q_ = epDict.lookupOrDefault<scalar>("q", 1e-2);
    }

    if (pcDict.found("density"))
    {
        const dictionary& densDict = pcDict.subDict("density");
        rhoRef_ = densDict.lookupOrDefault<scalar>("rhoRef", 1967.0);
        rhoSolid_ = densDict.lookupOrDefault<scalar>("rhoSolid", 1967.0);
        rhoLiquid_ = densDict.lookupOrDefault<scalar>("rhoLiquid", 1850.0);
    }
}

void Foam::pcmEnthalpyPorosityModel::correct()
{
    const volScalarField& T = thermo_.T();
    tmp<volScalarField> tCp = thermo_.Cp();
    tmp<volScalarField> tKappa = thermo_.kappa();
    const volScalarField& CpField = tCp();
    const volScalarField& KappaField = tKappa();

    forAll(T, cellI)
    {
        scalar Tcell = T[cellI];
        scalar alphaL = 0.0;

        if (Tcell <= Tlm_)
        {
            alphaL = 0.0;
            phaseState_[cellI] = 0.0;
        }
        else if (Tcell >= Tum_)
        {
            alphaL = 1.0;
            phaseState_[cellI] = 2.0;
        }
        else
        {
            alphaL = (Tcell - Tlm_) / (Tum_ - Tlm_);
            alphaL = max(0.0, min(1.0, alphaL));
            phaseState_[cellI] = 1.0;
        }

        liquidFraction_[cellI] = alphaL;
        rho_[cellI] = (1.0 - alphaL) * rhoSolid_ + alphaL * rhoLiquid_;
        Cp_[cellI] = CpField[cellI];
        k_[cellI] = KappaField[cellI];
    }

    liquidFraction_.correctBoundaryConditions();
    phaseState_.correctBoundaryConditions();
    rho_.correctBoundaryConditions();
    Cp_.correctBoundaryConditions();
    k_.correctBoundaryConditions();
}

Foam::tmp<Foam::volScalarField> Foam::pcmEnthalpyPorosityModel::latentHeatSource() const
{
    tmp<volScalarField> tSource
    (
        volScalarField::New
        (
            "latentHeatSource",
            mesh_,
            dimensionedScalar("zero", dimEnergy/dimTime/dimVolume, 0.0)
        )
    );

    volScalarField& source = tSource.ref();
    tmp<volScalarField> tdAlphaLdt = fvc::ddt(liquidFraction_);

    source = -rho_ * Lm_ * tdAlphaLdt();
    return tSource;
}

Foam::tmp<Foam::volVectorField> Foam::pcmEnthalpyPorosityModel::momentumSource() const
{
    tmp<volVectorField> tSource
    (
        volVectorField::New
        (
            "momentumSource",
            mesh_,
            dimensionedVector("zero", dimForce/dimVolume, Zero)
        )
    );

    if (mesh_.foundObject<volVectorField>("U"))
    {
        const volVectorField& U = mesh_.lookupObject<volVectorField>("U");
        volVectorField& source = tSource.ref();

        forAll(liquidFraction_, cellI)
        {
            scalar alphaL = max(q_, liquidFraction_[cellI]);
            scalar damping = -Cu_ * sqr(1.0 - alphaL) / (pow3(alphaL) + q_);
            source[cellI] = damping * U[cellI];
        }
    }

    return tSource;
}

// ************************************************************************* //
