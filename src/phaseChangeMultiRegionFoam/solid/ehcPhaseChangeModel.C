/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Website:  www.openfoam.com
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

#include "ehcPhaseChangeModel.H"
#include "addToRunTimeSelectionTable.H"
#include "calculatedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(ehcPhaseChangeModel, 0);
    addToRunTimeSelectionTable(phaseChangeModel, ehcPhaseChangeModel, dictionary);
}

// * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ehcPhaseChangeModel::ehcPhaseChangeModel
(
    const fvMesh& mesh,
    const solidThermo& thermo
)
:
    phaseChangeModel(mesh, thermo),
    Tlm_(303.15),
    Tum_(313.15),
    Lm_(163000.0),
    Tlf_(303.15),
    Tuf_(313.15),
    Lf_(163000.0),
    hysteresisActive_(true),
    reversalTol_(1e-6),
    densityModel_("linear"),
    rhoRef_(1967.0),
    rhoSolid_(1967.0),
    rhoLiquid_(1850.0),
    thermoMode_("thermo"),
    Cps_(1980.0),
    Cpl_(2320.0),
    ks_(0.50),
    kl_(0.47),
    CpEff_
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
    heatingTrajectory_
    (
        IOobject("heatingTrajectory", mesh.time().timeName(), mesh, IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
        mesh,
        dimensionedScalar("heatingTrajectory", dimless, 1.0),
        "calculated"
    ),
    phaseFraction_old_
    (
        IOobject("phaseFraction_old", mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
        phaseFraction_
    ),
    heatingTrajectory_old_
    (
        IOobject("heatingTrajectory_old", mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
        heatingTrajectory_
    )
{
    active_ = true;
    readDict();

    // Fresh start initialization for phaseFraction_ if not read from disk
    if (!phaseFractionRestored_)
    {
        const volScalarField& Tinit = thermo_.T();
        forAll(phaseFraction_, cellI)
        {
            scalar Tval = Tinit[cellI];
            scalar alpha0 = (Tval <= Tlm_) ? 0.0 : ((Tval >= Tum_) ? 1.0 : (Tval - Tlm_) / (Tum_ - Tlm_));
            phaseFraction_[cellI] = alpha0;
        }
    }

    // Restart handling: if phaseFraction was read but heatingTrajectory was missing, infer trajectory
    if (phaseFractionRestored_ && !heatingTrajectory_.headerOk())
    {
        forAll(heatingTrajectory_, cellI)
        {
            heatingTrajectory_[cellI] = (phaseFraction_[cellI] > 0.5) ? 1.0 : 0.0;
        }
        forAll(heatingTrajectory_.boundaryFieldRef(), patchi)
        {
            fvPatchScalarField& pHeat = heatingTrajectory_.boundaryFieldRef()[patchi];
            const fvPatchScalarField& pAlpha = phaseFraction_.boundaryField()[patchi];
            forAll(pHeat, facei)
            {
                pHeat[facei] = (pAlpha[facei] > 0.5) ? 1.0 : 0.0;
            }
        }
    }

    // Clamp initial phaseFraction_ internal and boundary patch values to [0, 1]
    forAll(phaseFraction_, cellI)
    {
        phaseFraction_[cellI] = max(0.0, min(1.0, phaseFraction_[cellI]));
    }
    forAll(phaseFraction_.boundaryFieldRef(), patchi)
    {
        fvPatchScalarField& pAlpha = phaseFraction_.boundaryFieldRef()[patchi];
        forAll(pAlpha, facei)
        {
            pAlpha[facei] = max(0.0, min(1.0, pAlpha[facei]));
        }
    }

    // Clamp heatingTrajectory_ internal and boundary patch values to 0.0 or 1.0
    forAll(heatingTrajectory_, cellI)
    {
        heatingTrajectory_[cellI] = (heatingTrajectory_[cellI] > 0.5) ? 1.0 : 0.0;
    }
    forAll(heatingTrajectory_.boundaryFieldRef(), patchi)
    {
        fvPatchScalarField& pHeat = heatingTrajectory_.boundaryFieldRef()[patchi];
        forAll(pHeat, facei)
        {
            pHeat[facei] = (pHeat[facei] > 0.5) ? 1.0 : 0.0;
        }
    }

    phaseFraction_old_ = phaseFraction_;
    heatingTrajectory_old_ = heatingTrajectory_;

    updateHistory();
    correct();

    // Validate that all non-constraint patches are 'calculated' to prevent
    // restart files with zeroGradient/fixedValue from silently overwriting
    // values computed in the boundary face loop.
    auto checkCalculatedPatches = [&](const volScalarField& f)
    {
        forAll(f.boundaryField(), pI)
        {
            const fvPatchScalarField& pf = f.boundaryField()[pI];
            if
            (
                pf.type() != calculatedFvPatchScalarField::typeName
             && !polyPatch::constraintType(pf.patch().patch().type())
            )
            {
                FatalErrorInFunction
                    << "Field '" << f.name() << "' patch '" << pf.patch().name()
                    << "' has type '" << pf.type() << "' but must be 'calculated'.\n"
                    << "Remove or correct the patch entry in the restart file."
                    << exit(FatalError);
            }
        }
    };
    checkCalculatedPatches(phaseFraction_);
    checkCalculatedPatches(heatingTrajectory_);
    checkCalculatedPatches(phaseState_);
    checkCalculatedPatches(CpEff_);
    checkCalculatedPatches(rho_);
    checkCalculatedPatches(k_);
}


// * * * * * * * * * * * * * * Private Functions * * * * * * * * * * * * * * //

void Foam::ehcPhaseChangeModel::readDict()
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

    readConvectionDict(pcDict);

    // Default freezing bounds inherit from melting bounds unless overridden
    Tlf_ = Tlm_;
    Tuf_ = Tum_;
    Lf_ = Lm_;

    if (pcDict.found("freezing"))
    {
        const dictionary& freezeDict = pcDict.subDict("freezing");
        Tlf_ = freezeDict.lookupOrDefault<scalar>("T_lowerBound", Tlf_);
        Tuf_ = freezeDict.lookupOrDefault<scalar>("T_upperBound", Tuf_);
        Lf_ = freezeDict.lookupOrDefault<scalar>("latentHeat", Lf_);
    }

    // Input Parameter Validation
    if (Tum_ <= Tlm_)
    {
        FatalIOErrorInFunction(phaseChangeDict)
            << "Melting upper bound T_upperBound (" << Tum_
            << " K) must be > T_lowerBound (" << Tlm_ << " K)"
            << exit(FatalIOError);
    }
    if (Tuf_ <= Tlf_)
    {
        FatalIOErrorInFunction(phaseChangeDict)
            << "Freezing upper bound T_upperBound (" << Tuf_
            << " K) must be > T_lowerBound (" << Tlf_ << " K)"
            << exit(FatalIOError);
    }
    if (Lm_ < 0.0 || Lf_ < 0.0)
    {
        FatalIOErrorInFunction(phaseChangeDict)
            << "Latent heat (Lm=" << Lm_ << ", Lf=" << Lf_ << ") must be >= 0"
            << exit(FatalIOError);
    }

    if (mag(Lm_ - Lf_) > 1e-6)
    {
        WarningInFunction
            << "Unequal melting latent heat (" << Lm_
            << " J/kg) and freezing latent heat (" << Lf_
            << " J/kg) specified for region " << mesh_.name()
            << ". Closed heating-cooling cycles will produce net enthalpy hysteresis."
            << endl;
    }

    if (pcDict.found("hysteresis"))
    {
        const dictionary& hysDict = pcDict.subDict("hysteresis");
        hysteresisActive_ = hysDict.lookupOrDefault<bool>("active", true);
        reversalTol_ = hysDict.lookupOrDefault<scalar>("reversalTol", 1e-6);
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

    tmp<volScalarField> tCpThermoInit = thermo_.Cp();
    const volScalarField& CpFieldInit = tCpThermoInit();
    scalar CpThermoInit = CpFieldInit.primitiveField().size() > 0 ? CpFieldInit.primitiveField()[0] : 1000.0;

    tmp<volScalarField> tKappaThermoInit = thermo_.kappa();
    const volScalarField& kFieldInit = tKappaThermoInit();
    scalar kThermoInit = kFieldInit.primitiveField().size() > 0 ? kFieldInit.primitiveField()[0] : 1.0;

    thermoMode_ = "thermo";
    Cps_ = CpThermoInit;
    Cpl_ = CpThermoInit;
    ks_ = kThermoInit;
    kl_ = kThermoInit;

    if (pcDict.found("thermophysical"))
    {
        const dictionary& thermoDict = pcDict.subDict("thermophysical");
        thermoMode_ = thermoDict.lookupOrDefault<word>("mode", "thermo");

        Cps_ = thermoDict.lookupOrDefault<scalar>("CpSolid", CpThermoInit);
        Cpl_ = thermoDict.lookupOrDefault<scalar>("CpLiquid", CpThermoInit);
        ks_ = thermoDict.lookupOrDefault<scalar>("kSolid", kThermoInit);
        kl_ = thermoDict.lookupOrDefault<scalar>("kLiquid", kThermoInit);
    }

    if (Cps_ <= 0.0 || Cpl_ <= 0.0)
    {
        FatalIOErrorInFunction(phaseChangeDict)
            << "Heat capacity (CpSolid=" << Cps_ << ", CpLiquid=" << Cpl_ << ") must be > 0"
            << exit(FatalIOError);
    }
    if (ks_ <= 0.0 || kl_ <= 0.0)
    {
        FatalIOErrorInFunction(phaseChangeDict)
            << "Thermal conductivity (kSolid=" << ks_ << ", kLiquid=" << kl_ << ") must be > 0"
            << exit(FatalIOError);
    }
    if (densityModel_ != "thermo" && densityModel_ != "constant" && densityModel_ != "linear")
    {
        FatalIOErrorInFunction(phaseChangeDict)
            << "Unsupported density model '" << densityModel_ << "'. Expected 'thermo', 'constant', or 'linear'."
            << exit(FatalIOError);
    }
    if (thermoMode_ != "thermo" && thermoMode_ != "custom")
    {
        FatalIOErrorInFunction(phaseChangeDict)
            << "Unsupported thermophysical mode '" << thermoMode_ << "'. Expected 'thermo' or 'custom'."
            << exit(FatalIOError);
    }

    if (hysteresisActive_ && mag(Cps_ - Cpl_) > 1e-6)
    {
        Info<< "    EHC Info: Unequal Cp (Cps=" << Cps_ << ", Cpl=" << Cpl_
            << " J/(kg.K)) under stateful hysteresis produces path-dependent sensible enthalpy loops across closed thermal cycles." << endl;
    }

    Info<< "    EHC Phase Change Parameters loaded for region " << mesh_.name() << ":" << nl
        << "      Melting Range: [" << Tlm_ << " - " << Tum_ << "] K, Latent Heat: " << Lm_ << " J/kg" << nl
        << "      Freezing Range: [" << Tlf_ << " - " << Tuf_ << "] K, Latent Heat: " << Lf_ << " J/kg" << nl
        << "      Hysteresis Active: " << (hysteresisActive_ ? "true" : "false") << nl
        << "      Density Model: " << densityModel_ << " (rhoS=" << rhoSolid_ << ", rhoL=" << rhoLiquid_ << ")" << endl;
}

Foam::scalar Foam::ehcPhaseChangeModel::integralAlphaMelt(scalar T, scalar alpha_old) const
{
    scalar Tstar = Tlm_ + alpha_old * (Tum_ - Tlm_);
    if (T <= Tstar)
    {
        return alpha_old * T;
    }
    else if (T < Tum_)
    {
        scalar I_Tstar = alpha_old * Tstar;
        scalar I_seg = ((T - Tlm_) * (T - Tlm_) - (Tstar - Tlm_) * (Tstar - Tlm_)) / (2.0 * (Tum_ - Tlm_));
        return I_Tstar + I_seg;
    }
    else
    {
        scalar I_Tstar = alpha_old * Tstar;
        scalar I_mush = ((Tum_ - Tlm_) * (Tum_ - Tlm_) - (Tstar - Tlm_) * (Tstar - Tlm_)) / (2.0 * (Tum_ - Tlm_));
        return I_Tstar + I_mush + (T - Tum_);
    }
}

Foam::scalar Foam::ehcPhaseChangeModel::integralAlphaFreeze(scalar T, scalar alpha_old) const
{
    scalar Tstar = Tlf_ + alpha_old * (Tuf_ - Tlf_);
    if (T <= Tlf_)
    {
        return 0.0;
    }
    else if (T < Tstar)
    {
        return (T - Tlf_) * (T - Tlf_) / (2.0 * (Tuf_ - Tlf_));
    }
    else
    {
        scalar I_mush = (Tstar - Tlf_) * (Tstar - Tlf_) / (2.0 * (Tuf_ - Tlf_));
        return I_mush + alpha_old * (T - Tstar);
    }
}


// * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * * //

void Foam::ehcPhaseChangeModel::updateHistory()
{
    phaseFraction_old_ = phaseFraction_;
    heatingTrajectory_old_ = heatingTrajectory_;
}


void Foam::ehcPhaseChangeModel::correct()
{
    const volScalarField& T = thermo_.T();
    const volScalarField& Told = thermo_.T().oldTime();
    tmp<volScalarField> tCpThermo = thermo_.Cp();
    tmp<volScalarField> tKappaThermo = thermo_.kappa();
    tmp<volScalarField> tRhoThermo = thermo_.rho();
    const volScalarField& CpField = tCpThermo();
    const volScalarField& KappaField = tKappaThermo();
    const volScalarField& rhoField = tRhoThermo();

    forAll(T, cellI)
    {
        scalar Tcell = T[cellI];
        scalar ToldCell = Told[cellI];
        scalar deltaTStep = Tcell - ToldCell;

        bool isHeating = (heatingTrajectory_old_[cellI] > 0.5);

        if (!hysteresisActive_)
        {
            isHeating = true;
            heatingTrajectory_[cellI] = 1.0;
        }
        else
        {
            if (deltaTStep > reversalTol_)
            {
                isHeating = true;
                heatingTrajectory_[cellI] = 1.0;
            }
            else if (deltaTStep < -reversalTol_)
            {
                isHeating = false;
                heatingTrajectory_[cellI] = 0.0;
            }
            else
            {
                heatingTrajectory_[cellI] = heatingTrajectory_old_[cellI];
            }
        }

        scalar L = isHeating ? Lm_ : Lf_;

        scalar alpha_m_curr = (Tcell <= Tlm_) ? 0.0 : ((Tcell >= Tum_) ? 1.0 : (Tcell - Tlm_) / (Tum_ - Tlm_));
        scalar alpha_f_curr = (Tcell <= Tlf_) ? 0.0 : ((Tcell >= Tuf_) ? 1.0 : (Tcell - Tlf_) / (Tuf_ - Tlf_));

        scalar alphaL_prev = phaseFraction_old_[cellI];
        if (!hysteresisActive_)
        {
            alphaL_prev = (ToldCell <= Tlm_) ? 0.0 : ((ToldCell >= Tum_) ? 1.0 : (ToldCell - Tlm_) / (Tum_ - Tlm_));
        }

        scalar alphaL_curr = 0.0;
        if (!hysteresisActive_)
        {
            alphaL_curr = alpha_m_curr;
        }
        else if (isHeating)
        {
            alphaL_curr = max(alphaL_prev, alpha_m_curr);
        }
        else
        {
            alphaL_curr = min(alphaL_prev, alpha_f_curr);
        }

        alphaL_curr = max(0.0, min(1.0, alphaL_curr));

        if (alphaL_curr <= 0.0) phaseState_[cellI] = 0.0;
        else if (alphaL_curr >= 1.0) phaseState_[cellI] = 2.0;
        else phaseState_[cellI] = isHeating ? 1.0 : 3.0;

        phaseFraction_[cellI] = alphaL_curr;

        scalar cpsVal = Cps_;
        scalar cplVal = Cpl_;
        scalar ksVal = ks_;
        scalar klVal = kl_;

        if (thermoMode_ == "thermo")
        {
            cpsVal = CpField[cellI];
            cplVal = CpField[cellI];
            ksVal = KappaField[cellI];
            klVal = KappaField[cellI];
        }

        scalar alpha_old_eval = hysteresisActive_ ? alphaL_prev : 0.0;
        scalar avgAlpha = 0.5 * (alphaL_curr + alphaL_prev);
        if (mag(deltaTStep) > reversalTol_)
        {
            if (!hysteresisActive_ || isHeating)
            {
                avgAlpha = (integralAlphaMelt(Tcell, alpha_old_eval) - integralAlphaMelt(ToldCell, alpha_old_eval)) / deltaTStep;
            }
            else
            {
                avgAlpha = (integralAlphaFreeze(Tcell, alpha_old_eval) - integralAlphaFreeze(ToldCell, alpha_old_eval)) / deltaTStep;
            }
        }
        avgAlpha = max(0.0, min(1.0, avgAlpha));

        scalar cpBase = (1.0 - avgAlpha) * cpsVal + avgAlpha * cplVal;

        scalar deltaCp = 0.0;
        if (mag(deltaTStep) > reversalTol_)
        {
            deltaCp = L * (alphaL_curr - alphaL_prev) / deltaTStep;
        }
        else
        {
            if (isHeating && Tcell > Tlm_ && Tcell < Tum_ && alphaL_prev <= alpha_m_curr + reversalTol_)
            {
                deltaCp = Lm_ / (Tum_ - Tlm_);
            }
            else if (!isHeating && Tcell > Tlf_ && Tcell < Tuf_ && alphaL_prev >= alpha_f_curr - reversalTol_)
            {
                deltaCp = Lf_ / (Tuf_ - Tlf_);
            }
        }

        CpEff_[cellI] = cpBase + max(0.0, deltaCp);

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
            rho_[cellI] = (1.0 - alphaL_curr) * rhoSolid_ + alphaL_curr * rhoLiquid_;
        }

        k_[cellI] = (1.0 - alphaL_curr) * ksVal + alphaL_curr * klVal;
    }

    // Update boundary patch face values
    forAll(T.boundaryField(), patchi)
    {
        const fvPatchScalarField& pT = T.boundaryField()[patchi];
        const fvPatchScalarField& pTold = Told.boundaryField()[patchi];
        fvPatchScalarField& pHeatTraj = heatingTrajectory_.boundaryFieldRef()[patchi];
        const fvPatchScalarField& pHeatTrajOld = heatingTrajectory_old_.boundaryField()[patchi];
        fvPatchScalarField& pPhaseFraction = phaseFraction_.boundaryFieldRef()[patchi];
        const fvPatchScalarField& pPhaseFractionOld = phaseFraction_old_.boundaryField()[patchi];
        fvPatchScalarField& pPhaseState = phaseState_.boundaryFieldRef()[patchi];
        fvPatchScalarField& pCpEff = CpEff_.boundaryFieldRef()[patchi];
        fvPatchScalarField& pRho = rho_.boundaryFieldRef()[patchi];
        fvPatchScalarField& pK = k_.boundaryFieldRef()[patchi];

        const fvPatchScalarField& pCpThermo = CpField.boundaryField()[patchi];
        const fvPatchScalarField& pKappaThermo = KappaField.boundaryField()[patchi];
        const fvPatchScalarField& pRhoThermo = rhoField.boundaryField()[patchi];

        forAll(pT, facei)
        {
            scalar Tface = pT[facei];
            scalar ToldFace = pTold[facei];
            scalar deltaTStep = Tface - ToldFace;

            bool isHeating = (pHeatTrajOld[facei] > 0.5);

            if (!hysteresisActive_)
            {
                isHeating = true;
                pHeatTraj[facei] = 1.0;
            }
            else
            {
                if (deltaTStep > reversalTol_)
                {
                    isHeating = true;
                    pHeatTraj[facei] = 1.0;
                }
                else if (deltaTStep < -reversalTol_)
                {
                    isHeating = false;
                    pHeatTraj[facei] = 0.0;
                }
                else
                {
                    pHeatTraj[facei] = pHeatTrajOld[facei];
                }
            }

            scalar L = isHeating ? Lm_ : Lf_;

            scalar alpha_m_curr = (Tface <= Tlm_) ? 0.0 : ((Tface >= Tum_) ? 1.0 : (Tface - Tlm_) / (Tum_ - Tlm_));
            scalar alpha_f_curr = (Tface <= Tlf_) ? 0.0 : ((Tface >= Tuf_) ? 1.0 : (Tface - Tlf_) / (Tuf_ - Tlf_));

            scalar alphaL_prev = pPhaseFractionOld[facei];
            if (!hysteresisActive_)
            {
                alphaL_prev = (ToldFace <= Tlm_) ? 0.0 : ((ToldFace >= Tum_) ? 1.0 : (ToldFace - Tlm_) / (Tum_ - Tlm_));
            }

            scalar alphaL_curr = 0.0;
            if (!hysteresisActive_)
            {
                alphaL_curr = alpha_m_curr;
            }
            else if (isHeating)
            {
                alphaL_curr = max(alphaL_prev, alpha_m_curr);
            }
            else
            {
                alphaL_curr = min(alphaL_prev, alpha_f_curr);
            }

            alphaL_curr = max(0.0, min(1.0, alphaL_curr));

            if (alphaL_curr <= 0.0) pPhaseState[facei] = 0.0;
            else if (alphaL_curr >= 1.0) pPhaseState[facei] = 2.0;
            else pPhaseState[facei] = isHeating ? 1.0 : 3.0;

            pPhaseFraction[facei] = alphaL_curr;

            scalar cpsVal = Cps_;
            scalar cplVal = Cpl_;
            scalar ksVal = ks_;
            scalar klVal = kl_;

            if (thermoMode_ == "thermo")
            {
                cpsVal = pCpThermo[facei];
                cplVal = pCpThermo[facei];
                ksVal = pKappaThermo[facei];
                klVal = pKappaThermo[facei];
            }

            scalar alpha_old_eval = hysteresisActive_ ? alphaL_prev : 0.0;
            scalar avgAlpha = 0.5 * (alphaL_curr + alphaL_prev);
            if (mag(deltaTStep) > reversalTol_)
            {
                if (!hysteresisActive_ || isHeating)
                {
                    avgAlpha = (integralAlphaMelt(Tface, alpha_old_eval) - integralAlphaMelt(ToldFace, alpha_old_eval)) / deltaTStep;
                }
                else
                {
                    avgAlpha = (integralAlphaFreeze(Tface, alpha_old_eval) - integralAlphaFreeze(ToldFace, alpha_old_eval)) / deltaTStep;
                }
            }
            avgAlpha = max(0.0, min(1.0, avgAlpha));

            scalar cpBase = (1.0 - avgAlpha) * cpsVal + avgAlpha * cplVal;

            scalar deltaCp = 0.0;
            if (mag(deltaTStep) > reversalTol_)
            {
                deltaCp = L * (alphaL_curr - alphaL_prev) / deltaTStep;
            }
            else
            {
                if (isHeating && Tface > Tlm_ && Tface < Tum_ && alphaL_prev <= alpha_m_curr + reversalTol_)
                {
                    deltaCp = Lm_ / (Tum_ - Tlm_);
                }
                else if (!isHeating && Tface > Tlf_ && Tface < Tuf_ && alphaL_prev >= alpha_f_curr - reversalTol_)
                {
                    deltaCp = Lf_ / (Tuf_ - Tlf_);
                }
            }

            pCpEff[facei] = cpBase + max(0.0, deltaCp);

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
                pRho[facei] = (1.0 - alphaL_curr) * rhoSolid_ + alphaL_curr * rhoLiquid_;
            }

            pK[facei] = (1.0 - alphaL_curr) * ksVal + alphaL_curr * klVal;
        }
    }

    CpEff_.correctBoundaryConditions();
    rho_.correctBoundaryConditions();
    k_.correctBoundaryConditions();
}

// ************************************************************************* //
