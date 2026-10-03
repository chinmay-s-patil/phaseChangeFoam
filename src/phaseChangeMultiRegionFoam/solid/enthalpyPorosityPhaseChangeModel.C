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
    rhoRef_(0.0),
    rhoSolid_(0.0),
    rhoLiquid_(0.0),
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

    if (!pcDict.found("melting"))
    {
        FatalIOErrorInFunction(phaseChangeDict)
            << "Missing required sub-dictionary 'melting' in 'phaseChange' for region "
            << mesh_.name()
            << exit(FatalIOError);
    }
    const dictionary& meltDict = pcDict.subDict("melting");
    if (!meltDict.found("T_lowerBound") || !meltDict.found("T_upperBound") || !meltDict.found("latentHeat"))
    {
        FatalIOErrorInFunction(phaseChangeDict)
            << "Sub-dictionary 'melting' in 'phaseChange' for region " << mesh_.name()
            << " must contain 'T_lowerBound', 'T_upperBound', and 'latentHeat'."
            << exit(FatalIOError);
    }
    Tlm_ = meltDict.get<scalar>("T_lowerBound");
    Tum_ = meltDict.get<scalar>("T_upperBound");
    Lm_ = meltDict.get<scalar>("latentHeat");

    if (Tum_ <= Tlm_)
    {
        FatalIOErrorInFunction(phaseChangeDict)
            << "Melting upper bound T_upperBound (" << Tum_
            << " K) must be > T_lowerBound (" << Tlm_ << " K)"
            << exit(FatalIOError);
    }
    if (Lm_ < 0.0)
    {
        FatalIOErrorInFunction(phaseChangeDict)
            << "Latent heat Lm (" << Lm_ << ") must be >= 0"
            << exit(FatalIOError);
    }

    if (pcDict.found("enthalpyPorosity"))
    {
        const dictionary& epDict = pcDict.subDict("enthalpyPorosity");
        beta_ = epDict.lookupOrDefault<scalar>("beta", 50e-6);
        Cu_ = epDict.lookupOrDefault<scalar>("Cu", 1e5);
        q_ = epDict.lookupOrDefault<scalar>("q", 1e-2);
    }

    tmp<volScalarField> tRhoThermoInit = thermo_.rho();
    const volScalarField& rhoFieldInit = tRhoThermoInit();
    scalar rhoThermoInit = rhoFieldInit.primitiveField().size() > 0 ? rhoFieldInit.primitiveField()[0] : 1000.0;
    densityModel_ = "thermo";
    rhoRef_ = rhoThermoInit;
    rhoSolid_ = rhoThermoInit;
    rhoLiquid_ = rhoThermoInit;

    if (pcDict.found("density"))
    {
        const dictionary& densDict = pcDict.subDict("density");
        densityModel_ = densDict.lookupOrDefault<word>("model", "thermo");
        rhoRef_ = densDict.lookupOrDefault<scalar>("rhoRef", rhoThermoInit);
        rhoSolid_ = densDict.lookupOrDefault<scalar>("rhoSolid", rhoThermoInit);
        rhoLiquid_ = densDict.lookupOrDefault<scalar>("rhoLiquid", rhoThermoInit);
    }

    if (densityModel_ != "thermo" && densityModel_ != "constant" && densityModel_ != "linear")
    {
        FatalIOErrorInFunction(phaseChangeDict)
            << "Unsupported density model '" << densityModel_ << "'. Expected 'thermo', 'constant', or 'linear'."
            << exit(FatalIOError);
    }
}

void Foam::enthalpyPorosityPhaseChangeModel::correct()
{
    const volScalarField& T = thermo_.T();
    tmp<volScalarField> tCp = thermo_.Cp();
    tmp<volScalarField> tKappa = thermo_.kappa();
    tmp<volScalarField> tRhoThermo = thermo_.rho();
    const volScalarField& CpField = tCp();
    const volScalarField& KappaField = tKappa();
    const volScalarField& rhoField = tRhoThermo();

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
        if (densityModel_ == "thermo")
        {
            rho_[cellI] = rhoField[cellI];
        }
        else if (densityModel_ == "constant")
        {
            rho_[cellI] = rhoRef_;
        }
        else
        {
            rho_[cellI] = (1.0 - alphaL) * rhoSolid_ + alphaL * rhoLiquid_;
        }
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
        const fvPatchScalarField& pRhoThermo = rhoField.boundaryField()[patchi];

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
            if (densityModel_ == "thermo")
            {
                pRho[facei] = pRhoThermo[facei];
            }
            else if (densityModel_ == "constant")
            {
                pRho[facei] = rhoRef_;
            }
            else
            {
                pRho[facei] = (1.0 - alphaL) * rhoSolid_ + alphaL * rhoLiquid_;
            }
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
    const volScalarField& T = thermo_.T();
    const volScalarField& h = thermo_.he();
    const scalar dAlpha = 1.0/(Tum_ - Tlm_);

    forAll(T, cellI)
    {
        scalar Tcell = T[cellI];
        scalar alphaOld = phaseFraction_old_[cellI];
        scalar rhoL_dt = rho_[cellI] * Lm_ * rDeltaT;

        if (Tcell <= Tlm_)
        {
            source[cellI] = rhoL_dt * alphaOld;
        }
        else if (Tcell >= Tum_)
        {
            source[cellI] = - rhoL_dt * (1.0 - alphaOld);
        }
        else
        {
            scalar SpVal = rhoL_dt * dAlpha / max(Cp_[cellI], 1e-10);
            scalar h_lm = h[cellI] - Cp_[cellI] * (Tcell - Tlm_);
            source[cellI] = SpVal * h_lm - rhoL_dt * alphaOld;
        }
    }

    return tSource;
}

Foam::tmp<Foam::volScalarField> Foam::enthalpyPorosityPhaseChangeModel::latentHeatSp() const
{
    return tmp<volScalarField>::New
    (
        IOobject("latentHeatSp", mesh_.time().timeName(), mesh_),
        mesh_,
        dimensionedScalar("zero", dimMass/dimVolume/dimTime, 0.0)
    );
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
