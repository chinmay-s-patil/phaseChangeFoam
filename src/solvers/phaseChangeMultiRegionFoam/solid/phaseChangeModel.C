/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Website:  www.openfoam.com
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

#include "phaseChangeModel.H"
#include "calculatedFvPatchFields.H"
#include <cctype>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(phaseChangeModel, 0);
    defineRunTimeSelectionTable(phaseChangeModel, dictionary);

    const Enum<phaseChangeModel::direction> phaseChangeModel::directionNames
    ({
        { direction::both, "both" },
        { direction::forward, "forward" },
        { direction::reverse, "reverse" }
    });
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseChangeModel::phaseChangeModel
(
    const fvMesh& mesh,
    const basicThermo& thermo,
    bool suppressConvection
)
:
    mesh_(mesh),
    thermo_(thermo),
    active_(false),
    dir_(direction::both),
    Tlm_(300.0),
    Tum_(310.0),
    Lm_(1.0e5),
    Tlf_(300.0),
    Tuf_(310.0),
    Lf_(1.0e5),
    hysteresisActive_(false),
    reversalTol_(1e-6),
    densityModel_("thermo"),
    rhoRef_(1000.0),
    rhoSolid_(1000.0),
    rhoLiquid_(1000.0),
    thermoMode_("thermo"),
    Cps_(1000.0),
    Cpl_(1000.0),
    ks_(1.0),
    kl_(1.0),
    phaseFraction_
    (
        IOobject
        (
            "phaseFraction",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("phaseFraction", dimless, 0.0),
        "calculated"
    ),
    phaseState_
    (
        IOobject
        (
            "phaseState",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("phaseState", dimless, 0.0),
        "calculated"
    ),
    Cp_
    (
        IOobject
        (
            "CpEff",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        thermo.Cp(),
        "calculated"
    ),
    rho_
    (
        IOobject
        (
            "rhoEff",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        thermo.rho(),
        "calculated"
    ),
    k_
    (
        IOobject
        (
            "kEff",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        thermo.kappa(),
        "calculated"
    ),
    Su_
    (
        IOobject
        (
            "latentHeatSource",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimEnergy/dimTime/dimVolume, 0.0)
    ),
    Sp_
    (
        IOobject
        (
            "latentHeatSp",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimMass/dimVolume/dimTime, 0.0)
    ),
    heatingTrajectory_
    (
        IOobject
        (
            "heatingTrajectory",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("heatingTrajectory", dimless, 1.0),
        "calculated"
    ),
    T_reversal_
    (
        IOobject
        (
            "T_reversal",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("T_reversal", dimTemperature, 0.0),
        "calculated"
    ),
    phaseFraction_old_
    (
        IOobject
        (
            "phaseFraction_old",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        phaseFraction_
    ),
    heatingTrajectory_old_
    (
        IOobject
        (
            "heatingTrajectory_old",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        heatingTrajectory_
    ),
    T_reversal_old_
    (
        IOobject
        (
            "T_reversal_old",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("T_reversal_old", dimTemperature, 0.0),
        "calculated"
    ),
    phaseFractionRestored_(false),
    suppressConvection_(suppressConvection)
{
    phaseFractionRestored_ = phaseFraction_.headerOk();

    if (!T_reversal_.headerOk())
    {
        T_reversal_ == thermo.T();
    }

    // Backward compatibility for disk reads: if liquidFraction exists on disk but phaseFraction doesn't
    IOobject liquidFractionIO
    (
        "liquidFraction",
        mesh.time().timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    );

    if (!phaseFraction_.headerOk() && liquidFractionIO.typeHeaderOk<volScalarField>(true))
    {
        volScalarField lFrac(liquidFractionIO, mesh);
        phaseFraction_ == lFrac;
        phaseFractionRestored_ = true;
    }

    phaseFraction_old_ == phaseFraction_;
    heatingTrajectory_old_ == heatingTrajectory_;
    T_reversal_old_ == T_reversal_;

    readDict();
}


namespace
{
    void checkAllowedKeys
    (
        const dictionary& dict,
        const wordHashSet& allowedKeys
    )
    {
        for (const word& key : dict.toc())
        {
            if (!allowedKeys.found(key))
            {
                FatalIOErrorInFunction(dict)
                    << "Unknown key '" << key << "' in dictionary '" << dict.dictName()
                    << "'. Allowed keys: " << allowedKeys
                    << exit(FatalIOError);
            }
        }
    }
}


void Foam::phaseChangeModel::readDict()
{
    IOobject dictIO
    (
        "phaseChangeDict",
        mesh_.time().constant(),
        mesh_,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    if (!dictIO.typeHeaderOk<dictionary>(true))
    {
        return;
    }

    IOdictionary phaseChangeDict(dictIO);

    if (!phaseChangeDict.found("phaseChange"))
    {
        return;
    }

    const dictionary& pcDict = phaseChangeDict.subDict("phaseChange");
    checkAllowedKeys(pcDict, {"active", "phaseChangeMode", "type", "direction", "melting", "freezing", "forward", "reverse", "hysteresis", "density", "thermophysical", "thermo", "convection", "porosity", "buoyancy", "liquid", "vapor", "C_evap", "C_cond", "latentHeat", "Tsat", "enableTsatP", "pRef", "speciesName", "massSource", "energySource", "speciesSource"});

    word modeName = "none";
    if (pcDict.found("phaseChangeMode"))
    {
        modeName = pcDict.get<word>("phaseChangeMode");
    }
    else if (pcDict.found("type"))
    {
        modeName = pcDict.get<word>("type");
    }

    auto toLowerStr = [](const word& w)
    {
        std::string s(w);
        for (char& c : s)
        {
            c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
        }
        return s;
    };

    word modeLower = toLowerStr(modeName);
    if (modeLower != "ehc" && modeLower != "ehcdarcy" && modeLower != "enthalpyporosity")
    {
        active_ = false;
        return;
    }

    active_ = pcDict.lookupOrDefault<bool>("active", true);

    if (!active_)
    {
        return;
    }

    dir_ = directionNames.getOrDefault("direction", pcDict, direction::both);

    if (pcDict.found("forward") && pcDict.found("melting"))
    {
        FatalIOErrorInFunction(pcDict)
            << "Cannot specify both 'forward' and deprecated 'melting' sub-dictionaries in phaseChangeDict for region "
            << mesh_.name() << exit(FatalIOError);
    }
    if (pcDict.found("reverse") && pcDict.found("freezing"))
    {
        FatalIOErrorInFunction(pcDict)
            << "Cannot specify both 'reverse' and deprecated 'freezing' sub-dictionaries in phaseChangeDict for region "
            << mesh_.name() << exit(FatalIOError);
    }

    word forwardKey = pcDict.found("forward") ? "forward" : (pcDict.found("melting") ? "melting" : "");
    if (forwardKey == "melting")
    {
        WarningInFunction
            << "Key 'melting' in phaseChange dictionary is deprecated. Use 'forward' instead." << endl;
    }

    word reverseKey = pcDict.found("reverse") ? "reverse" : (pcDict.found("freezing") ? "freezing" : "");
    if (reverseKey == "freezing")
    {
        WarningInFunction
            << "Key 'freezing' in phaseChange dictionary is deprecated. Use 'reverse' instead." << endl;
    }

    if (dir_ == direction::reverse)
    {
        // Reverse / Freezing sub-dictionary is mandatory in reverse direction mode
        if (reverseKey.empty())
        {
            FatalIOErrorInFunction(pcDict)
                << "Mandatory 'reverse' (or legacy 'freezing') block missing in phaseChangeDict for region "
                << mesh_.name() << " with direction reverse." << exit(FatalIOError);
        }

        if (!forwardKey.empty())
        {
            WarningInFunction
                << "Key '" << forwardKey << "' in phaseChange dictionary is ignored when direction is 'reverse'." << endl;
        }

        const dictionary& freezeDict = pcDict.subDict(reverseKey);
        checkAllowedKeys(freezeDict, {"T_lowerBound", "T_upperBound", "latentHeat"});
        Tlf_ = freezeDict.get<scalar>("T_lowerBound");
        Tuf_ = freezeDict.get<scalar>("T_upperBound");
        Lf_  = freezeDict.get<scalar>("latentHeat");

        if (Tuf_ <= Tlf_)
        {
            FatalIOErrorInFunction(freezeDict)
                << "T_upperBound (" << Tuf_ << ") must be strictly greater than "
                << "T_lowerBound (" << Tlf_ << ") in reverse/freezing block for region "
                << mesh_.name() << exit(FatalIOError);
        }

        // Forward defaults to reverse values
        Tlm_ = Tlf_;
        Tum_ = Tuf_;
        Lm_  = Lf_;

        if (!forwardKey.empty())
        {
            const dictionary& meltDict = pcDict.subDict(forwardKey);
            checkAllowedKeys(meltDict, {"T_lowerBound", "T_upperBound", "latentHeat"});
            Tlm_ = meltDict.lookupOrDefault<scalar>("T_lowerBound", Tlf_);
            Tum_ = meltDict.lookupOrDefault<scalar>("T_upperBound", Tuf_);
            Lm_  = meltDict.lookupOrDefault<scalar>("latentHeat", Lf_);

            if (Tum_ <= Tlm_)
            {
                FatalIOErrorInFunction(meltDict)
                    << "T_upperBound (" << Tum_ << ") must be strictly greater than "
                    << "T_lowerBound (" << Tlm_ << ") in forward/melting block for region "
                    << mesh_.name() << exit(FatalIOError);
            }
        }
    }
    else
    {
        // Forward / Melting sub-dictionary is mandatory in both or forward direction mode
        if (forwardKey.empty())
        {
            FatalIOErrorInFunction(pcDict)
                << "Mandatory 'forward' (or legacy 'melting') block missing in phaseChangeDict for region "
                << mesh_.name() << exit(FatalIOError);
        }

        if (dir_ == direction::forward && !reverseKey.empty())
        {
            WarningInFunction
                << "Key '" << reverseKey << "' in phaseChange dictionary is ignored when direction is 'forward'." << endl;
        }

        const dictionary& meltDict = pcDict.subDict(forwardKey);
        checkAllowedKeys(meltDict, {"T_lowerBound", "T_upperBound", "latentHeat"});
        Tlm_ = meltDict.get<scalar>("T_lowerBound");
        Tum_ = meltDict.get<scalar>("T_upperBound");
        Lm_  = meltDict.get<scalar>("latentHeat");

        if (Tum_ <= Tlm_)
        {
            FatalIOErrorInFunction(meltDict)
                << "T_upperBound (" << Tum_ << ") must be strictly greater than "
                << "T_lowerBound (" << Tlm_ << ") in forward/melting block for region "
                << mesh_.name() << exit(FatalIOError);
        }

        // Reverse defaults to forward values
        Tlf_ = Tlm_;
        Tuf_ = Tum_;
        Lf_  = Lm_;

        if (!reverseKey.empty())
        {
            const dictionary& freezeDict = pcDict.subDict(reverseKey);
            checkAllowedKeys(freezeDict, {"T_lowerBound", "T_upperBound", "latentHeat"});
            Tlf_ = freezeDict.lookupOrDefault<scalar>("T_lowerBound", Tlm_);
            Tuf_ = freezeDict.lookupOrDefault<scalar>("T_upperBound", Tum_);
            Lf_  = freezeDict.lookupOrDefault<scalar>("latentHeat", Lm_);

            if (Tuf_ <= Tlf_)
            {
                FatalIOErrorInFunction(freezeDict)
                    << "T_upperBound (" << Tuf_ << ") must be strictly greater than "
                    << "T_lowerBound (" << Tlf_ << ") in reverse/freezing block for region "
                    << mesh_.name() << exit(FatalIOError);
            }
        }
    }

    // Hysteresis sub-dictionary
    if (pcDict.found("hysteresis"))
    {
        const dictionary& hysDict = pcDict.subDict("hysteresis");
        checkAllowedKeys(hysDict, {"active", "reversalTolerance", "reversalTol"});
        hysteresisActive_ = hysDict.lookupOrDefault<bool>("active", false);
        if (hysDict.found("reversalTol") && !hysDict.found("reversalTolerance"))
        {
            WarningInFunction
                << "Key 'reversalTol' in hysteresis dictionary is deprecated. Use 'reversalTolerance' instead." << endl;
            reversalTol_ = hysDict.get<scalar>("reversalTol");
        }
        else
        {
            reversalTol_ = hysDict.lookupOrDefault<scalar>("reversalTolerance", 1e-6);
        }
    }

    // Density model
    if (pcDict.found("density"))
    {
        const dictionary& rhoDict = pcDict.subDict("density");
        checkAllowedKeys(rhoDict, {"model", "rhoRef", "rhoSolid", "rhoLiquid", "solid", "liquid", "allowNonConservativeDensity"});
        densityModel_ = rhoDict.lookupOrDefault<word>("model", "thermo");
        rhoRef_   = rhoDict.lookupOrDefault<scalar>("rhoRef", 1000.0);
        rhoSolid_ = rhoDict.getOrDefault<scalar>("rhoSolid", rhoDict.lookupOrDefault<scalar>("solid", rhoRef_));
        rhoLiquid_= rhoDict.getOrDefault<scalar>("rhoLiquid", rhoDict.lookupOrDefault<scalar>("liquid", rhoRef_));
        allowNonConservativeDensity_ = rhoDict.lookupOrDefault<bool>("allowNonConservativeDensity", false);

        if (densityModel_ == "linear" && mag(rhoSolid_ - rhoLiquid_) > 1e-4 && !allowNonConservativeDensity_)
        {
            FatalIOErrorInFunction(pcDict)
                << "Density model 'linear' with rhoSolid (" << rhoSolid_
                << ") != rhoLiquid (" << rhoLiquid_
                << ") in a fixed-volume solid domain is non-conservative for energy storage.\n"
                << "Set 'allowNonConservativeDensity true;' in phaseChange.density dictionary to proceed, "
                << "or set density model to 'thermo' (constant density)."
                << exit(FatalIOError);
        }
    }

    // Thermophysical model (custom vs thermo)
    if (pcDict.found("thermophysical") || pcDict.found("thermo"))
    {
        const dictionary& thDict = pcDict.found("thermophysical")
            ? pcDict.subDict("thermophysical")
            : pcDict.subDict("thermo");
        checkAllowedKeys(thDict, {"mode", "CpSolid", "CpLiquid", "kSolid", "kLiquid", "Cps", "Cpl", "ks", "kl"});
        thermoMode_ = thDict.lookupOrDefault<word>("mode", "thermo");
        if (thermoMode_ == "custom")
        {
            Cps_ = thDict.getOrDefault<scalar>("CpSolid", thDict.lookupOrDefault<scalar>("Cps", 1000.0));
            Cpl_ = thDict.getOrDefault<scalar>("CpLiquid", thDict.lookupOrDefault<scalar>("Cpl", 1000.0));
            ks_  = thDict.getOrDefault<scalar>("kSolid", thDict.lookupOrDefault<scalar>("ks", 1.0));
            kl_  = thDict.getOrDefault<scalar>("kLiquid", thDict.lookupOrDefault<scalar>("kl", 1.0));

            scalar maxCpDiff = gMax(mag(thermo_.Cp()().primitiveField() - Cps_));
            if (maxCpDiff > 1e-4)
            {
                FatalIOErrorInFunction(pcDict)
                    << "In custom thermo mode, thermophysicalProperties Cp must match CpSolid (" << Cps_
                    << ") across all cells (max diff = " << maxCpDiff << ")."
                    << exit(FatalIOError);
            }
        }
    }

    readConvectionDict(pcDict);
}


void Foam::phaseChangeModel::readConvectionDict(const dictionary& pcDict)
{
    if (pcDict.found("convection"))
    {
        const dictionary& convDict = pcDict.subDict("convection");
        checkAllowedKeys(convDict, {"suppress"});
        suppressConvection_ =
            convDict.getOrDefault<bool>("suppress", suppressConvection_);
    }

    if (!suppressConvection_ && !mesh_.foundObject<volVectorField>("U"))
    {
        FatalIOErrorInFunction(pcDict)
            << "phaseChange.convection.suppress = false requested for region "
            << mesh_.name() << ", but this region does not solve velocity field U. "
            << "Set suppress true or use EHCDarcy mode on a fluid region."
            << exit(FatalIOError);
    }
}


Foam::phaseChangeModel::StateResult Foam::phaseChangeModel::evalState
(
    scalar Tc,
    scalar Told,
    scalar aOld,
    scalar trajOld,
    scalar TrevOld
) const
{
    StateResult res;

    // Trajectory tracking with last reversal temperature tracking
    scalar traj = trajOld;
    scalar Trev = TrevOld;

    if (traj > 0.5) // Heating branch
    {
        if (Tc > Trev)
        {
            Trev = Tc;
        }
        else if (Trev - Tc > reversalTol_)
        {
            traj = 0.0; // Cooling reversal
            Trev = Tc;
        }
    }
    else // Cooling branch
    {
        if (Tc < Trev)
        {
            Trev = Tc;
        }
        else if (Tc - Trev > reversalTol_)
        {
            traj = 1.0; // Heating reversal
            Trev = Tc;
        }
    }

    // Direction flag enforcement (override trajectory if forward or reverse)
    if (dir_ == direction::forward)
    {
        traj = 1.0;
    }
    else if (dir_ == direction::reverse)
    {
        traj = 0.0;
    }

    res.traj = traj;
    res.Trev = Trev;

    // Phase fraction & state machine calculation
    scalar aVal = 0.0;
    bool isPlateau = false;

    scalar a_melt = (Tum_ > Tlm_) ? clamp((Tc - Tlm_) / (Tum_ - Tlm_), 0.0, 1.0) : (Tc >= Tum_ ? 1.0 : 0.0);
    scalar a_freeze = (Tuf_ > Tlf_) ? clamp((Tc - Tlf_) / (Tuf_ - Tlf_), 0.0, 1.0) : (Tc >= Tuf_ ? 1.0 : 0.0);

    if (traj > 0.5) // Heating branch
    {
        if (hysteresisActive_ && a_melt < aOld)
        {
            aVal = aOld;
            isPlateau = true;
        }
        else
        {
            aVal = a_melt;
        }
    }
    else // Cooling branch
    {
        if (hysteresisActive_ && a_freeze > aOld)
        {
            aVal = aOld;
            isPlateau = true;
        }
        else
        {
            aVal = a_freeze;
        }
    }

    // Direction restriction clamp: forward (irreversible melting/evaporation) or reverse (irreversible freezing/condensation)
    bool clamped = false;
    if (dir_ == direction::forward && aVal < aOld)
    {
        aVal = aOld;
        clamped = true;
    }
    else if (dir_ == direction::reverse && aVal > aOld)
    {
        aVal = aOld;
        clamped = true;
    }

    isPlateau = isPlateau || clamped;

    // Compute phaseState AFTER direction clamp and hysteresis plateau check
    scalar stateVal = 0.0; // 0 = solid, 1 = melting, 2 = liquid, 3 = freezing
    if (aVal <= 0.0) stateVal = 0.0;
    else if (aVal >= 1.0) stateVal = 2.0;
    else stateVal = (traj > 0.5) ? 1.0 : 3.0;

    res.aVal = aVal;
    res.stateVal = stateVal;
    res.isPlateau = isPlateau;

    // Determine dAlpha_dT for Newton tangent stabilization
    scalar dAlpha_dT = 0.0;
    if (!isPlateau)
    {
        if (traj > 0.5) // Heating branch
        {
            if ((Tc > Tlm_ && Tc < Tum_) || (Tc >= Tum_ && aOld < 1.0))
            {
                dAlpha_dT = (Tum_ > Tlm_) ? (1.0 / (Tum_ - Tlm_)) : 0.0;
            }
        }
        else // Cooling branch
        {
            if ((Tc > Tlf_ && Tc < Tuf_) || (Tc <= Tlf_ && aOld > 0.0))
            {
                dAlpha_dT = (Tuf_ > Tlf_) ? (1.0 / (Tuf_ - Tlf_)) : 0.0;
            }
        }
    }
    res.dAlpha_dT = dAlpha_dT;

    return res;
}


Foam::scalar Foam::phaseChangeModel::hSens(scalar T, scalar traj) const
{
    scalar Tl = (traj > 0.5) ? Tlm_ : Tlf_;
    scalar Tu = (traj > 0.5) ? Tum_ : Tuf_;

    if (T <= Tl)
    {
        return Cps_ * T;
    }
    else if (T >= Tu)
    {
        scalar h_l = Cps_ * Tl;
        scalar h_mush = (0.5 * (Cps_ + Cpl_)) * (Tu - Tl);
        return h_l + h_mush + Cpl_ * (T - Tu);
    }
    else
    {
        scalar h_l = Cps_ * Tl;
        scalar dT = T - Tl;
        scalar dTum = Tu - Tl;
        if (dTum > 1e-8)
        {
            return h_l + Cps_ * dT + 0.5 * (Cpl_ - Cps_) * sqr(dT) / dTum;
        }
        else
        {
            return h_l;
        }
    }
}


Foam::scalar Foam::phaseChangeModel::deltaHSens
(
    scalar Tnew,
    scalar Told,
    scalar aOld,
    scalar trajNew,
    scalar trajOld
) const
{
    scalar dT = Tnew - Told;
    if (mag(dT) < 1e-12)
    {
        return 0.0;
    }

    if (hysteresisActive_ || dir_ != direction::both)
    {
        if (trajNew > 0.5) // Heating branch
        {
            scalar Texit = Tlm_ + aOld * (Tum_ - Tlm_);
            if (Texit < Tlm_) Texit = Tlm_;
            if (Texit > Tum_) Texit = Tum_;

            if (Told < Texit)
            {
                if (Tnew <= Texit)
                {
                    scalar Cp_plateau = (1.0 - aOld) * Cps_ + aOld * Cpl_;
                    return Cp_plateau * dT;
                }
                else
                {
                    scalar Cp_plateau = (1.0 - aOld) * Cps_ + aOld * Cpl_;
                    scalar dH1 = Cp_plateau * (Texit - Told);
                    scalar dH2 = hSens(Tnew, 1.0) - hSens(Texit, 1.0);
                    return dH1 + dH2;
                }
            }
            else // Told >= Texit (Cooling step on heating branch, e.g. direction::forward)
            {
                if (Tnew >= Texit)
                {
                    return hSens(Tnew, 1.0) - hSens(Told, 1.0);
                }
                else
                {
                    scalar Cp_plateau = (1.0 - aOld) * Cps_ + aOld * Cpl_;
                    scalar dH1 = hSens(Texit, 1.0) - hSens(Told, 1.0);
                    scalar dH2 = Cp_plateau * (Tnew - Texit);
                    return dH1 + dH2;
                }
            }
        }
        else // Cooling branch
        {
            scalar Texit = Tlf_ + aOld * (Tuf_ - Tlf_);
            if (Texit < Tlf_) Texit = Tlf_;
            if (Texit > Tuf_) Texit = Tuf_;

            if (Told > Texit)
            {
                if (Tnew >= Texit)
                {
                    scalar Cp_plateau = (1.0 - aOld) * Cps_ + aOld * Cpl_;
                    return Cp_plateau * dT;
                }
                else
                {
                    scalar Cp_plateau = (1.0 - aOld) * Cps_ + aOld * Cpl_;
                    scalar dH1 = Cp_plateau * (Texit - Told);
                    scalar dH2 = hSens(Tnew, 0.0) - hSens(Texit, 0.0);
                    return dH1 + dH2;
                }
            }
            else // Told <= Texit (Heating step on cooling branch, e.g. direction::reverse)
            {
                if (Tnew <= Texit)
                {
                    return hSens(Tnew, 0.0) - hSens(Told, 0.0);
                }
                else
                {
                    scalar Cp_plateau = (1.0 - aOld) * Cps_ + aOld * Cpl_;
                    scalar dH1 = hSens(Texit, 0.0) - hSens(Told, 0.0);
                    scalar dH2 = Cp_plateau * (Tnew - Texit);
                    return dH1 + dH2;
                }
            }
        }
    }

    return hSens(Tnew, trajNew) - hSens(Told, trajNew);
}


void Foam::phaseChangeModel::correct()
{
    if (!active_) return;

    const_cast<basicThermo&>(thermo_).correct();

    const volScalarField& T = thermo_.T();
    const volScalarField& T_old = thermo_.T().oldTime();
    scalar dt = mesh_.time().deltaTValue();
    if (dt <= 0.0) dt = 1.0;

    tmp<volScalarField> tCpThermo = thermo_.Cp();
    const scalarField& CpThermoCells = tCpThermo().primitiveField();
    tmp<volScalarField> tRhoThermo = thermo_.rho();
    const scalarField& rhoThermoCells = tRhoThermo().primitiveField();
    tmp<volScalarField> tKappaThermo = thermo_.kappa();
    const scalarField& kappaThermoCells = tKappaThermo().primitiveField();

    scalarField& alphaCells = phaseFraction_.primitiveFieldRef();
    scalarField& stateCells = phaseState_.primitiveFieldRef();
    scalarField& trajCells = heatingTrajectory_.primitiveFieldRef();
    scalarField& TrevCells = T_reversal_.primitiveFieldRef();
    scalarField& CpCells = Cp_.primitiveFieldRef();
    scalarField& rhoCells = rho_.primitiveFieldRef();
    scalarField& kCells = k_.primitiveFieldRef();
    scalarField& SuCells = Su_.primitiveFieldRef();
    scalarField& SpCells = Sp_.primitiveFieldRef();

    const scalarField& Tcells = T.primitiveField();
    const scalarField& ToldCells = T_old.primitiveField();
    const scalarField& alphaOldCells = phaseFraction_old_.primitiveField();
    const scalarField& trajOldCells = heatingTrajectory_old_.primitiveField();
    const scalarField& TrevOldCells = T_reversal_old_.primitiveField();

    forAll(Tcells, celli)
    {
        scalar Tc = Tcells[celli];
        scalar Told = ToldCells[celli];
        scalar aOld = alphaOldCells[celli];

        StateResult st = evalState
        (
            Tc,
            Told,
            aOld,
            trajOldCells[celli],
            TrevOldCells[celli]
        );

        trajCells[celli] = st.traj;
        TrevCells[celli] = st.Trev;
        alphaCells[celli] = st.aVal;
        stateCells[celli] = st.stateVal;

        // Sensible heat capacity Cp calculation (path-secant C_p,eff)
        scalar cpVal = CpThermoCells[celli];
        if (thermoMode_ == "custom")
        {
            scalar dT_step = Tc - Told;
            if (mag(dT_step) > 1e-8)
            {
                cpVal = deltaHSens(Tc, Told, aOld, st.traj, trajOldCells[celli]) / dT_step;
            }
            else
            {
                cpVal = (1.0 - st.aVal) * Cps_ + st.aVal * Cpl_;
            }
            kCells[celli] = (1.0 - st.aVal) * ks_ + st.aVal * kl_;
        }
        else
        {
            kCells[celli] = kappaThermoCells[celli];
        }
        CpCells[celli] = cpVal;

        // Density calculation (path-integrated rho_bar for energy terms)
        scalar rhoVal = rhoThermoCells[celli];
        scalar rhoBar = rhoThermoCells[celli];
        if (densityModel_ == "linear")
        {
            rhoVal = (1.0 - st.aVal) * rhoSolid_ + st.aVal * rhoLiquid_;
            scalar a_bar = 0.5 * (st.aVal + aOld);
            rhoBar = (1.0 - a_bar) * rhoSolid_ + a_bar * rhoLiquid_;
        }
        rhoCells[celli] = rhoVal;

        // Newton tangent slope Sp = + (rho_bar * L * dAlpha_dT) / (dt * Cp_thermo)
        // Explicit latent heat source Su = + (rho_bar * L * (aVal - aOld)) / dt
        // In energy equation hEqn: hEqn += Su + fvm::Sp(Sp, h) - Sp * h_k (on LHS), so Su > 0 acts as a sink during melting
        scalar L = (st.traj > 0.5) ? Lm_ : Lf_;
        scalar SpVal = 0.0;
        scalar cpThVal = CpThermoCells[celli];
        if (st.dAlpha_dT > 0.0 && cpThVal > 0.0)
        {
            SpVal = (rhoBar * L * st.dAlpha_dT) / (dt * cpThVal);
        }

        scalar SuVal = (rhoBar * L * (st.aVal - aOld)) / dt;

        SpCells[celli] = SpVal;
        SuCells[celli] = SuVal;
    }


    // Per-face boundary evaluations for non-constraint patches
    forAll(mesh_.boundaryMesh(), patchI)
    {
        const polyPatch& pp = mesh_.boundaryMesh()[patchI];
        if (polyPatch::constraintType(pp.type()))
        {
            continue;
        }

        const scalarField& Tface = T.boundaryField()[patchI];
        const scalarField& ToldFace = T_old.boundaryField()[patchI];
        const scalarField& alphaOldFace = phaseFraction_old_.boundaryField()[patchI];
        const scalarField& trajOldFace = heatingTrajectory_old_.boundaryField()[patchI];
        const scalarField& TrevOldFace = T_reversal_old_.boundaryField()[patchI];

        scalarField& alphaFace = phaseFraction_.boundaryFieldRef()[patchI];
        scalarField& stateFace = phaseState_.boundaryFieldRef()[patchI];
        scalarField& trajFace = heatingTrajectory_.boundaryFieldRef()[patchI];
        scalarField& TrevFace = T_reversal_.boundaryFieldRef()[patchI];
        scalarField& CpFace = Cp_.boundaryFieldRef()[patchI];
        scalarField& rhoFace = rho_.boundaryFieldRef()[patchI];
        scalarField& kFace = k_.boundaryFieldRef()[patchI];
        scalarField& SuFace = Su_.boundaryFieldRef()[patchI];
        scalarField& SpFace = Sp_.boundaryFieldRef()[patchI];

        const scalarField& CpThermoFace = tCpThermo().boundaryField()[patchI];
        const scalarField& rhoThermoFace = tRhoThermo().boundaryField()[patchI];
        const scalarField& kappaThermoFace = tKappaThermo().boundaryField()[patchI];

        forAll(Tface, faceI)
        {
            scalar Tf = Tface[faceI];
            scalar Toldf = ToldFace[faceI];
            scalar aOldf = alphaOldFace[faceI];

            StateResult stf = evalState
            (
                Tf,
                Toldf,
                aOldf,
                trajOldFace[faceI],
                TrevOldFace[faceI]
            );

            trajFace[faceI] = stf.traj;
            TrevFace[faceI] = stf.Trev;
            alphaFace[faceI] = stf.aVal;
            stateFace[faceI] = stf.stateVal;

            scalar cpValf = CpThermoFace[faceI];
            if (thermoMode_ == "custom")
            {
                scalar dT_stepf = Tf - Toldf;
                if (mag(dT_stepf) > 1e-8)
                {
                    cpValf = deltaHSens(Tf, Toldf, aOldf, stf.traj, trajOldFace[faceI]) / dT_stepf;
                }
                else
                {
                    cpValf = (1.0 - stf.aVal) * Cps_ + stf.aVal * Cpl_;
                }
                kFace[faceI] = (1.0 - stf.aVal) * ks_ + stf.aVal * kl_;
            }
            else
            {
                kFace[faceI] = kappaThermoFace[faceI];
            }
            CpFace[faceI] = cpValf;

            scalar rhoValf = rhoThermoFace[faceI];
            scalar rhoBarf = rhoThermoFace[faceI];
            if (densityModel_ == "linear")
            {
                rhoValf = (1.0 - stf.aVal) * rhoSolid_ + stf.aVal * rhoLiquid_;
                scalar a_barf = 0.5 * (stf.aVal + aOldf);
                rhoBarf = (1.0 - a_barf) * rhoSolid_ + a_barf * rhoLiquid_;
            }
            rhoFace[faceI] = rhoValf;

            scalar Lf_val = (stf.traj > 0.5) ? Lm_ : Lf_;
            scalar SpValf = 0.0;
            scalar cpThValf = CpThermoFace[faceI];
            if (stf.dAlpha_dT > 0.0 && cpThValf > 0.0)
            {
                SpValf = (rhoBarf * Lf_val * stf.dAlpha_dT) / (dt * cpThValf);
            }

            scalar SuValf = (rhoBarf * Lf_val * (stf.aVal - aOldf)) / dt;

            SpFace[faceI] = SpValf;
            SuFace[faceI] = SuValf;
        }
    }

    phaseFraction_.correctBoundaryConditions();
    phaseState_.correctBoundaryConditions();
    rho_.correctBoundaryConditions();
    Cp_.correctBoundaryConditions();
    k_.correctBoundaryConditions();
    Su_.correctBoundaryConditions();
    Sp_.correctBoundaryConditions();
    T_reversal_.correctBoundaryConditions();
}


void Foam::phaseChangeModel::updateHistory()
{
    if (active_)
    {
        phaseFraction_old_ == phaseFraction_;
        heatingTrajectory_old_ == heatingTrajectory_;
        T_reversal_old_ == T_reversal_;
    }
}


// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::phaseChangeModel> Foam::phaseChangeModel::New
(
    const fvMesh& mesh,
    const basicThermo& thermo
)
{
    IOobject dictIO
    (
        "phaseChangeDict",
        mesh.time().constant(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    );

    if (dictIO.typeHeaderOk<dictionary>(true))
    {
        IOdictionary phaseChangeDict(dictIO);

        bool active = phaseChangeDict.lookupOrDefault<bool>("active", true);

        if (active && phaseChangeDict.found("phaseChange"))
        {
            const dictionary& pcDict = phaseChangeDict.subDict("phaseChange");
            bool pcActive = pcDict.lookupOrDefault<bool>("active", true);

            if (pcActive)
            {
                word modeName = "none";
                if (pcDict.found("phaseChangeMode"))
                {
                    modeName = pcDict.get<word>("phaseChangeMode");
                }
                else if (pcDict.found("type"))
                {
                    word typeName = pcDict.get<word>("type");
                    auto cstrTest = dictionaryConstructorTablePtr_->cfind(typeName);
                    if (cstrTest.good())
                    {
                        modeName = typeName;
                    }
                    else
                    {
                        auto toLowerStr = [](const word& w)
                        {
                            std::string s(w);
                            for (char& c : s)
                            {
                                c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
                            }
                            return s;
                        };

                        for (auto iter = dictionaryConstructorTablePtr_->cbegin(); iter != dictionaryConstructorTablePtr_->cend(); ++iter)
                        {
                            if (toLowerStr(iter.key()) == toLowerStr(typeName))
                            {
                                modeName = iter.key();
                                break;
                            }
                        }
                    }
                }

                auto cstrIter = dictionaryConstructorTablePtr_->cfind(modeName);

                if (cstrIter.good() && modeName != "none")
                {
                    Info<< "    PCM Phase change ACTIVE for region " << mesh.name()
                        << " using mode: " << modeName << endl;
                    return autoPtr<phaseChangeModel>(cstrIter()(mesh, thermo));
                }
            }
        }
    }

    Info<< "    PCM Phase change INACTIVE for region " << mesh.name() << endl;

    auto cstrIter = dictionaryConstructorTablePtr_->cfind("none");
    if (cstrIter.good())
    {
        return autoPtr<phaseChangeModel>(cstrIter()(mesh, thermo));
    }

    return autoPtr<phaseChangeModel>(nullptr);
}

// ************************************************************************* //
