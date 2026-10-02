/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Website:  www.openfoam.com
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

#include "enthalpyPorosityPhaseChangeModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(enthalpyPorosityPhaseChangeModel, 0);
    addToRunTimeSelectionTable(phaseChangeModel, enthalpyPorosityPhaseChangeModel, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::enthalpyPorosityPhaseChangeModel::enthalpyPorosityPhaseChangeModel
(
    const fvMesh& mesh,
    const solidThermo& thermo
)
:
    phaseChangeModel(mesh, thermo),
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
        IOobject("CpEff", mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
        thermo.Cp()
    ),
    rho_
    (
        IOobject("rhoEff", mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
        thermo.rho()
    ),
    k_
    (
        IOobject("kEff", mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
        thermo.kappa()
    ),
    phaseFraction_old_
    (
        IOobject("phaseFraction_old", mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
        phaseFraction_
    )
{
    active_ = true;
    readDict();
    correct();
    phaseFraction_old_ = phaseFraction_;
}

void Foam::enthalpyPorosityPhaseChangeModel::readDict()
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

    readConvectionDict(pcDict);

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

void Foam::enthalpyPorosityPhaseChangeModel::correct()
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

        phaseFraction_[cellI] = alphaL;
        rho_[cellI] = (1.0 - alphaL) * rhoSolid_ + alphaL * rhoLiquid_;
        Cp_[cellI] = CpField[cellI];
        k_[cellI] = KappaField[cellI];
    }

    forAll(T.boundaryField(), patchi)
    {
        const fvPatchScalarField& pT = T.boundaryField()[patchi];
        fvPatchScalarField& pPhaseFraction = phaseFraction_.boundaryFieldRef()[patchi];
        fvPatchScalarField& pPhaseState = phaseState_.boundaryFieldRef()[patchi];
        fvPatchScalarField& pRho = rho_.boundaryFieldRef()[patchi];
        fvPatchScalarField& pCp = Cp_.boundaryFieldRef()[patchi];
        fvPatchScalarField& pK = k_.boundaryFieldRef()[patchi];

        const fvPatchScalarField& pCpThermo = CpField.boundaryField()[patchi];
        const fvPatchScalarField& pKappaThermo = KappaField.boundaryField()[patchi];

        forAll(pT, facei)
        {
            scalar Tface = pT[facei];
            scalar alphaL = 0.0;

            if (Tface <= Tlm_)
            {
                alphaL = 0.0;
                pPhaseState[facei] = 0.0;
            }
            else if (Tface >= Tum_)
            {
                alphaL = 1.0;
                pPhaseState[facei] = 2.0;
            }
            else
            {
                alphaL = (Tface - Tlm_) / (Tum_ - Tlm_);
                alphaL = max(0.0, min(1.0, alphaL));
                pPhaseState[facei] = 1.0;
            }

            pPhaseFraction[facei] = alphaL;
            pRho[facei] = (1.0 - alphaL) * rhoSolid_ + alphaL * rhoLiquid_;
            pCp[facei] = pCpThermo[facei];
            pK[facei] = pKappaThermo[facei];
        }
    }

    rho_.correctBoundaryConditions();
    Cp_.correctBoundaryConditions();
    k_.correctBoundaryConditions();
}

void Foam::enthalpyPorosityPhaseChangeModel::updateHistory()
{
    phaseFraction_old_ = phaseFraction_;
}


Foam::tmp<Foam::volScalarField> Foam::enthalpyPorosityPhaseChangeModel::latentHeatSource() const
{
    // Term is added to the LHS of the energy equation (hEqn += source), i.e.
    //   rho Cp dT/dt - div(k grad T) + rho L dalpha/dt = ...
    // so melting (dalpha/dt > 0) absorbs heat.
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
    const scalar rDeltaT = 1.0/mesh_.time().deltaTValue();

    source.primitiveFieldRef() =
        rho_.primitiveField()*Lm_
       *(phaseFraction_.primitiveField() - phaseFraction_old_.primitiveField())
       *rDeltaT;

    return tSource;
}

Foam::tmp<Foam::volVectorField> Foam::enthalpyPorosityPhaseChangeModel::momentumSource() const
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

        forAll(phaseFraction_, cellI)
        {
            scalar alphaL = max(q_, phaseFraction_[cellI]);
            scalar damping = -Cu_ * sqr(1.0 - alphaL) / (pow3(alphaL) + q_);
            source[cellI] = damping * U[cellI];
        }
    }

    return tSource;
}

// ************************************************************************* //
