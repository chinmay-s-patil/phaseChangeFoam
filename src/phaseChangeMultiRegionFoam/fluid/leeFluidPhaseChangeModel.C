/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
\*---------------------------------------------------------------------------*/

#include "leeFluidPhaseChangeModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(leeFluidPhaseChangeModel, 0);
    addToRunTimeSelectionTable
    (
        fluidPhaseChangeModel,
        leeFluidPhaseChangeModel,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::leeFluidPhaseChangeModel::leeFluidPhaseChangeModel
(
    const fvMesh& mesh,
    const rhoReactionThermo& thermo,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    fluidPhaseChangeModel(mesh, thermo, U, phi),
    liquidName_(dict_.getOrDefault<word>("liquid", "liquid")),
    vaporName_(dict_.getOrDefault<word>("vapor", "vapor")),
    liquidIndex_(thermo_.composition().species().find(liquidName_)),
    vaporIndex_(thermo_.composition().species().find(vaporName_)),
    timeIndex_(-1),
    C_evap_
    (
        "C_evap",
        dimless/dimTime,
        dict_.getOrDefault<scalar>("C_evap", 0.1)
    ),
    C_cond_
    (
        "C_cond",
        dimless/dimTime,
        dict_.getOrDefault<scalar>("C_cond", 0.1)
    ),
    latentHeat_
    (
        "latentHeat",
        dimEnergy/dimMass,
        dict_.getOrDefault<scalar>("latentHeat", 2.26e6)
    ),
    TsatRef_
    (
        "Tsat",
        dimTemperature,
        dict_.getOrDefault<scalar>("Tsat", 373.15)
    ),
    pRef_
    (
        "pRef",
        dimPressure,
        dict_.getOrDefault<scalar>("pRef", 101325)
    ),
    enableTsatP_(dict_.getOrDefault<bool>("enableTsatP", false)),
    mDot_
    (
        IOobject
        (
            "mDotLee",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimMass/dimVolume/dimTime, Zero)
    ),
    massSource_
    (
        IOobject
        (
            "massSourceLee",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimMass/dimVolume/dimTime, Zero)
    ),
    energySource_
    (
        IOobject
        (
            "energySourceLee",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimEnergy/dimVolume/dimTime, Zero)
    ),
    liquidSource_
    (
        IOobject
        (
            "liquidSourceLee",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimMass/dimVolume/dimTime, Zero)
    ),
    vaporSource_
    (
        IOobject
        (
            "vaporSourceLee",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimMass/dimVolume/dimTime, Zero)
    )
{
    active_ = true;

    if (active_)
    {
        if (liquidIndex_ == -1)
        {
            FatalIOErrorInFunction(dict_)
                << "liquid specie '" << liquidName_
                << "' not found in thermo composition. Available species are: "
                << thermo_.composition().species()
                << exit(FatalIOError);
        }
        if (vaporIndex_ == -1)
        {
            FatalIOErrorInFunction(dict_)
                << "vapor specie '" << vaporName_
                << "' not found in thermo composition. Available species are: "
                << thermo_.composition().species()
                << exit(FatalIOError);
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::leeFluidPhaseChangeModel::~leeFluidPhaseChangeModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::leeFluidPhaseChangeModel::massSource() const
{
    return massSource_;
}


Foam::tmp<Foam::volScalarField>
Foam::leeFluidPhaseChangeModel::speciesSource(const label specieIndex) const
{
    if (specieIndex == liquidIndex_ && liquidIndex_ != -1)
    {
        return liquidSource_;
    }
    else if (specieIndex == vaporIndex_ && vaporIndex_ != -1)
    {
        return vaporSource_;
    }
    else
    {
        return zeroSpeciesSource_;
    }
}


Foam::tmp<Foam::volScalarField>
Foam::leeFluidPhaseChangeModel::energySource() const
{
    return energySource_;
}


void Foam::leeFluidPhaseChangeModel::correct()
{
    if (!active_)
    {
        return;
    }

    const volScalarField& T_old = thermo_.T().oldTime();
    const volScalarField& p_old = thermo_.p().oldTime();
    
    const volScalarField& rho_old = mesh_.lookupObject<volScalarField>("rho").oldTime();

    // thermo_.Cp() returns a fresh tmp<volScalarField> without historic oldTime() storage.
    // Because correct() is called at oCorr == 0 before any equation solves (EEqn/YEqn),
    // thermo_.Cp() evaluates Cp directly at the old-time state (T_0, Y_0).
    tmp<volScalarField> tCp = thermo_.Cp();
    const volScalarField& CpField = tCp();

    const PtrList<volScalarField>& Y = thermo_.composition().Y();

    const scalar dt = mesh_.time().deltaTValue();
    const scalar Tsat0 = TsatRef_.value();
    const scalar L = max(latentHeat_.value(), scalar(1e-6));
    const scalar C_evap = C_evap_.value();
    const scalar C_cond = C_cond_.value();

    scalar R_vap = 461.5; // Default gas constant for water vapor (J/kg/K)
    if (vaporIndex_ != -1)
    {
        const scalar W_v = thermo_.composition().W(vaporIndex_);
        if (W_v > 1e-3)
        {
            R_vap = 8314.463 / W_v; // Universal gas constant 8314.463 J/(kmol K)
        }
    }

    forAll(mesh_.cells(), cellI)
    {
        scalar Tsat_c = Tsat0;
        if (enableTsatP_)
        {
            const scalar p_c = p_old[cellI];
            if (p_c > 1e-3 && pRef_.value() > 1e-3)
            {
                scalar invTsat = (1.0 / Tsat0) - (R_vap / L) * ::log(p_c / pRef_.value());
                if (invTsat > 1e-6)
                {
                    Tsat_c = 1.0 / invTsat;
                }
            }
        }

        const scalar T0 = T_old[cellI];
        const scalar rho0 = rho_old[cellI];
        const scalar Cp_c = max(CpField[cellI], scalar(1.0));
        const scalar Yl0 = (liquidIndex_ != -1) ? max(Y[liquidIndex_].oldTime()[cellI], scalar(0)) : scalar(0);
        const scalar Yv0 = (vaporIndex_ != -1) ? max(Y[vaporIndex_].oldTime()[cellI], scalar(0)) : scalar(0);

        scalar mDotVal = 0.0;
        scalar energySourceVal = 0.0;

        if (T0 > Tsat_c)
        {
            const scalar mDot_raw = C_evap * rho0 * Yl0 * max(T0 - Tsat_c, scalar(0)) / Tsat_c;
            const scalar massCap = (dt > 1e-12) ? (rho0 * Yl0 / dt) : mDot_raw;
            const scalar thermalCap = (dt > 1e-12) ? (rho0 * Cp_c * max(T0 - Tsat_c, scalar(0)) / (L * dt)) : mDot_raw;

            mDotVal = min(mDot_raw, min(massCap, thermalCap));
            energySourceVal = - mDotVal * L;
        }
        else if (T0 < Tsat_c)
        {
            const scalar mDot_raw = C_cond * rho0 * Yv0 * max(Tsat_c - T0, scalar(0)) / Tsat_c;
            const scalar massCap = (dt > 1e-12) ? (rho0 * Yv0 / dt) : mDot_raw;
            const scalar thermalCap = (dt > 1e-12) ? (rho0 * Cp_c * max(Tsat_c - T0, scalar(0)) / (L * dt)) : mDot_raw;

            mDotVal = - min(mDot_raw, min(massCap, thermalCap));
            energySourceVal = - mDotVal * L;
        }

        mDot_[cellI] = mDotVal;
        liquidSource_[cellI] = - mDotVal;
        vaporSource_[cellI] = mDotVal;
        energySource_[cellI] = energySourceVal;
    }
}


// ************************************************************************* //
