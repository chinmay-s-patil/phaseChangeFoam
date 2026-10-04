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
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseChangeModel::phaseChangeModel
(
    const fvMesh& mesh,
    const solidThermo& thermo,
    bool suppressConvection
)
:
    mesh_(mesh),
    thermo_(thermo),
    active_(false),
    Tlm_(300.0),
    Tum_(310.0),
    Lm_(1.0e5),
    Tlf_(300.0),
    Tuf_(310.0),
    Lf_(1.0e5),
    hysteresisActive_(false),
    reversalTol_(1e-4),
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
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("heatingTrajectory", dimless, 1.0),
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
    phaseFractionRestored_(false),
    suppressConvection_(suppressConvection)
{
    phaseFractionRestored_ = phaseFraction_.headerOk();

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

    readDict();
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
    active_ = pcDict.lookupOrDefault<bool>("active", true);

    if (!active_)
    {
        return;
    }

    // Melting sub-dictionary is mandatory
    if (!pcDict.found("melting"))
    {
        FatalIOErrorInFunction(pcDict)
            << "Mandatory 'melting' block missing in phaseChangeDict for region "
            << mesh_.name() << exit(FatalIOError);
    }
    const dictionary& meltDict = pcDict.subDict("melting");
    Tlm_ = meltDict.get<scalar>("T_lowerBound");
    Tum_ = meltDict.get<scalar>("T_upperBound");
    Lm_  = meltDict.get<scalar>("latentHeat");

    if (Tum_ <= Tlm_)
    {
        FatalIOErrorInFunction(meltDict)
            << "T_upperBound (" << Tum_ << ") must be strictly greater than "
            << "T_lowerBound (" << Tlm_ << ") in melting block for region "
            << mesh_.name() << exit(FatalIOError);
    }

    // Freezing sub-dictionary (defaults to melting values if hysteresis inactive)
    Tlf_ = Tlm_;
    Tuf_ = Tum_;
    Lf_  = Lm_;

    if (pcDict.found("freezing"))
    {
        const dictionary& freezeDict = pcDict.subDict("freezing");
        Tlf_ = freezeDict.lookupOrDefault<scalar>("T_lowerBound", Tlm_);
        Tuf_ = freezeDict.lookupOrDefault<scalar>("T_upperBound", Tum_);
        Lf_  = freezeDict.lookupOrDefault<scalar>("latentHeat", Lm_);
    }

    // Hysteresis sub-dictionary
    if (pcDict.found("hysteresis"))
    {
        const dictionary& hysDict = pcDict.subDict("hysteresis");
        hysteresisActive_ = hysDict.lookupOrDefault<bool>("active", false);
        reversalTol_ = hysDict.lookupOrDefault<scalar>("reversalTolerance", 1e-4);
    }

    // Density model
    if (pcDict.found("density"))
    {
        const dictionary& rhoDict = pcDict.subDict("density");
        densityModel_ = rhoDict.lookupOrDefault<word>("model", "thermo");
        rhoRef_   = rhoDict.lookupOrDefault<scalar>("rhoRef", 1000.0);
        rhoSolid_ = rhoDict.lookupOrDefault<scalar>("solid", rhoRef_);
        rhoLiquid_= rhoDict.lookupOrDefault<scalar>("liquid", rhoRef_);
    }

    // Thermophysical model (custom vs thermo)
    if (pcDict.found("thermophysical") || pcDict.found("thermo"))
    {
        const dictionary& thDict = pcDict.found("thermophysical")
            ? pcDict.subDict("thermophysical")
            : pcDict.subDict("thermo");
        thermoMode_ = thDict.lookupOrDefault<word>("mode", "thermo");
        if (thermoMode_ == "custom")
        {
            Cps_ = thDict.getOrDefault<scalar>("CpSolid", thDict.lookupOrDefault<scalar>("Cps", 1000.0));
            Cpl_ = thDict.getOrDefault<scalar>("CpLiquid", thDict.lookupOrDefault<scalar>("Cpl", 1000.0));
            ks_  = thDict.getOrDefault<scalar>("kSolid", thDict.lookupOrDefault<scalar>("ks", 1.0));
            kl_  = thDict.getOrDefault<scalar>("kLiquid", thDict.lookupOrDefault<scalar>("kl", 1.0));

            scalar thermoCpVal = thermo_.Cp()().primitiveField()[0];
            if (mag(thermoCpVal - Cps_) > 1e-4)
            {
                FatalIOErrorInFunction(pcDict)
                    << "In custom thermo mode, thermophysicalProperties Cp (" << thermoCpVal
                    << ") must match CpSolid (" << Cps_ << ") to ensure consistent boundary face enthalpy."
                    << exit(FatalIOError);
            }
        }
    }

    readConvectionDict(pcDict);
}


void Foam::phaseChangeModel::readConvectionDict(const dictionary& pcDict)
{
    bool defaultSuppress = suppressConvection_;
    if (pcDict.found("convection"))
    {
        suppressConvection_ =
            pcDict.subDict("convection").getOrDefault<bool>("suppress", suppressConvection_);
    }

    if (defaultSuppress && !suppressConvection_)
    {
        FatalIOErrorInFunction(pcDict)
            << "phaseChange.convection.suppress = false requested for region "
            << mesh_.name() << ", but this model is configured for solid/conduction-only "
            << "regions. Set suppress true or use enthalpyPorosity mode."
            << exit(FatalIOError);
    }
}


void Foam::phaseChangeModel::correct()
{
    if (!active_) return;

    const_cast<solidThermo&>(thermo_).correct();

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
    scalarField& CpCells = Cp_.primitiveFieldRef();
    scalarField& rhoCells = rho_.primitiveFieldRef();
    scalarField& kCells = k_.primitiveFieldRef();
    scalarField& SuCells = Su_.primitiveFieldRef();
    scalarField& SpCells = Sp_.primitiveFieldRef();

    const scalarField& Tcells = T.primitiveField();
    const scalarField& ToldCells = T_old.primitiveField();
    const scalarField& alphaOldCells = phaseFraction_old_.primitiveField();
    const scalarField& trajOldCells = heatingTrajectory_old_.primitiveField();

    forAll(Tcells, celli)
    {
        scalar Tc = Tcells[celli];
        scalar Told = ToldCells[celli];
        scalar aOld = alphaOldCells[celli];

        // Trajectory tracking
        scalar dT = Tc - Told;
        scalar traj = trajOldCells[celli];
        if (dT > reversalTol_)
        {
            traj = 1.0; // Heating
        }
        else if (dT < -reversalTol_)
        {
            traj = 0.0; // Cooling
        }
        trajCells[celli] = traj;

        // Phase fraction & state machine calculation
        scalar aVal = 0.0;
        scalar stateVal = 0.0; // 0 = solid, 1 = melting, 2 = liquid, 3 = freezing
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

            if (aVal <= 0.0) stateVal = 0.0;
            else if (aVal >= 1.0) stateVal = 2.0;
            else stateVal = 1.0;
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

            if (aVal <= 0.0) stateVal = 0.0;
            else if (aVal >= 1.0) stateVal = 2.0;
            else stateVal = 3.0;
        }

        alphaCells[celli] = aVal;
        stateCells[celli] = stateVal;

        // Determine dAlpha_dT for Newton tangent stabilization
        scalar dAlpha_dT = 0.0;
        if (!isPlateau)
        {
            if (aVal != aOld || (Tc >= min(Tlm_, Tlf_) && Tc <= max(Tum_, Tuf_)))
            {
                if (traj > 0.5)
                {
                    dAlpha_dT = (Tum_ > Tlm_) ? (1.0 / (Tum_ - Tlm_)) : 0.0;
                }
                else
                {
                    dAlpha_dT = (Tuf_ > Tlf_) ? (1.0 / (Tuf_ - Tlf_)) : 0.0;
                }
            }
        }

        // Sensible heat capacity Cp calculation (state-fraction C_p(alpha))
        scalar cpVal = CpThermoCells[celli];
        if (thermoMode_ == "custom")
        {
            cpVal = (1.0 - aVal) * Cps_ + aVal * Cpl_;
            kCells[celli] = (1.0 - aVal) * ks_ + aVal * kl_;
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
            rhoVal = (1.0 - aVal) * rhoSolid_ + aVal * rhoLiquid_;
            scalar a_bar = 0.5 * (aVal + aOld);
            rhoBar = (1.0 - a_bar) * rhoSolid_ + a_bar * rhoLiquid_;
        }
        rhoCells[celli] = rhoVal;

        // Newton tangent slope Sp = + (rho_bar * L * dAlpha_dT) / (dt * cp)
        // Physical latent heat source Su_phys = + (rho_bar * L * (aVal - aOld)) / dt
        scalar L = (traj > 0.5) ? Lm_ : Lf_;
        scalar SpVal = 0.0;
        if (dAlpha_dT > 0.0 && cpVal > 0.0)
        {
            SpVal = (rhoBar * L * dAlpha_dT) / (dt * cpVal);
        }

        scalar SuVal = (rhoBar * L * (aVal - aOld)) / dt;

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

        scalarField& alphaFace = phaseFraction_.boundaryFieldRef()[patchI];
        scalarField& stateFace = phaseState_.boundaryFieldRef()[patchI];
        scalarField& trajFace = heatingTrajectory_.boundaryFieldRef()[patchI];
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

            scalar dTf = Tf - Toldf;
            scalar trajf = trajOldFace[faceI];
            if (dTf > reversalTol_)
            {
                trajf = 1.0;
            }
            else if (dTf < -reversalTol_)
            {
                trajf = 0.0;
            }
            trajFace[faceI] = trajf;

            scalar aValf = 0.0;
            scalar stateValf = 0.0;
            bool isPlateauf = false;

            scalar a_meltf = (Tum_ > Tlm_) ? clamp((Tf - Tlm_) / (Tum_ - Tlm_), 0.0, 1.0) : (Tf >= Tum_ ? 1.0 : 0.0);
            scalar a_freezef = (Tuf_ > Tlf_) ? clamp((Tf - Tlf_) / (Tuf_ - Tlf_), 0.0, 1.0) : (Tf >= Tuf_ ? 1.0 : 0.0);

            if (trajf > 0.5) // Heating branch
            {
                if (hysteresisActive_ && a_meltf < aOldf)
                {
                    aValf = aOldf;
                    isPlateauf = true;
                }
                else
                {
                    aValf = a_meltf;
                }

                if (aValf <= 0.0) stateValf = 0.0;
                else if (aValf >= 1.0) stateValf = 2.0;
                else stateValf = 1.0;
            }
            else // Cooling branch
            {
                if (hysteresisActive_ && a_freezef > aOldf)
                {
                    aValf = aOldf;
                    isPlateauf = true;
                }
                else
                {
                    aValf = a_freezef;
                }

                if (aValf <= 0.0) stateValf = 0.0;
                else if (aValf >= 1.0) stateValf = 2.0;
                else stateValf = 3.0;
            }

            alphaFace[faceI] = aValf;
            stateFace[faceI] = stateValf;

            scalar dAlpha_dTf = 0.0;
            if (!isPlateauf)
            {
                if (aValf != aOldf || (Tf >= min(Tlm_, Tlf_) && Tf <= max(Tum_, Tuf_)))
                {
                    if (trajf > 0.5)
                    {
                        dAlpha_dTf = (Tum_ > Tlm_) ? (1.0 / (Tum_ - Tlm_)) : 0.0;
                    }
                    else
                    {
                        dAlpha_dTf = (Tuf_ > Tlf_) ? (1.0 / (Tuf_ - Tlf_)) : 0.0;
                    }
                }
            }

            scalar cpValf = CpThermoFace[faceI];
            if (thermoMode_ == "custom")
            {
                cpValf = (1.0 - aValf) * Cps_ + aValf * Cpl_;
                kFace[faceI] = (1.0 - aValf) * ks_ + aValf * kl_;
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
                rhoValf = (1.0 - aValf) * rhoSolid_ + aValf * rhoLiquid_;
                scalar a_barf = 0.5 * (aValf + aOldf);
                rhoBarf = (1.0 - a_barf) * rhoSolid_ + a_barf * rhoLiquid_;
            }
            rhoFace[faceI] = rhoValf;

            scalar Lf_val = (trajf > 0.5) ? Lm_ : Lf_;
            scalar SpValf = 0.0;
            if (dAlpha_dTf > 0.0 && cpValf > 0.0)
            {
                SpValf = (rhoBarf * Lf_val * dAlpha_dTf) / (dt * cpValf);
            }

            scalar SuValf = (rhoBarf * Lf_val * (aValf - aOldf)) / dt;

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
}


void Foam::phaseChangeModel::updateHistory()
{
    if (active_)
    {
        phaseFraction_old_ == phaseFraction_;
        heatingTrajectory_old_ == heatingTrajectory_;
    }
}


// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::phaseChangeModel> Foam::phaseChangeModel::New
(
    const fvMesh& mesh,
    const solidThermo& thermo
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
                word modeName = pcDict.lookupOrDefault<word>("type", pcDict.lookupOrDefault<word>("phaseChangeMode", "ehc"));

                Info<< "    Phase change ACTIVE for region " << mesh.name()
                    << " using mode: " << modeName << endl;

                auto cstrIter = dictionaryConstructorTablePtr_->cfind(modeName);

                if (!cstrIter.good())
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
                        if (toLowerStr(iter.key()) == toLowerStr(modeName))
                        {
                            cstrIter = iter;
                            break;
                        }
                    }
                }

                if (!cstrIter.good())
                {
                    FatalErrorInFunction
                        << "Unknown phaseChange mode " << modeName
                        << " for region " << mesh.name() << nl << nl
                        << "Valid options are: "
                        << dictionaryConstructorTablePtr_->toc()
                        << exit(FatalError);
                }

                return autoPtr<phaseChangeModel>(cstrIter()(mesh, thermo));
            }
        }
    }

    Info<< "    Phase change INACTIVE for region " << mesh.name() << endl;

    auto cstrIter = dictionaryConstructorTablePtr_->cfind("none");
    if (cstrIter.good())
    {
        return autoPtr<phaseChangeModel>(cstrIter()(mesh, thermo));
    }

    return autoPtr<phaseChangeModel>(nullptr);
}

// ************************************************************************* //
