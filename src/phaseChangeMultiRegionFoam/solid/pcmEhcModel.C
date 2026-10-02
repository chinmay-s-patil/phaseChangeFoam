/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Website:  www.openfoam.com
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

#include "pcmEhcModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(pcmEhcModel, 0);
    addToRunTimeSelectionTable(pcmPhaseChangeModel, pcmEhcModel, dictionary);
}

// * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pcmEhcModel::pcmEhcModel
(
    const fvMesh& mesh,
    const solidThermo& thermo
)
:
    pcmPhaseChangeModel(mesh, thermo),
    Tlm_(303.15),
    Tum_(313.15),
    Lm_(163000.0),
    Tlf_(303.15),
    Tuf_(313.15),
    Lf_(163000.0),
    hysteresisActive_(true),
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
        IOobject("rhoPCM", mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
        thermo.rho()
    ),
    k_
    (
        IOobject("kPCM", mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
        thermo.kappa()
    ),
    heatingTrajectory_
    (
        IOobject("heatingTrajectory", mesh.time().timeName(), mesh, IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
        mesh,
        dimensionedScalar("heatingTrajectory", dimless, 1.0)
    ),
    liquidFraction_old_
    (
        IOobject("liquidFraction_old", mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
        liquidFraction_
    ),
    heatingTrajectory_old_
    (
        IOobject("heatingTrajectory_old", mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
        heatingTrajectory_
    )
{
    active_ = true;
    readDict();

    // Fresh start initialization for liquidFraction_ if not read from disk
    if (!liquidFraction_.headerOk())
    {
        const volScalarField& Tinit = thermo_.T();
        forAll(liquidFraction_, cellI)
        {
            scalar Tval = Tinit[cellI];
            scalar alpha0 = (Tval <= Tlm_) ? 0.0 : ((Tval >= Tum_) ? 1.0 : (Tval - Tlm_) / (Tum_ - Tlm_));
            liquidFraction_[cellI] = alpha0;
        }
    }
    liquidFraction_.correctBoundaryConditions();
    liquidFraction_old_ = liquidFraction_;
    heatingTrajectory_old_ = heatingTrajectory_;

    updateHistory();
    correct();
}


// * * * * * * * * * * * * * * Private Functions * * * * * * * * * * * * * * //

void Foam::pcmEhcModel::readDict()
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

    if (pcDict.found("freezing"))
    {
        const dictionary& freezeDict = pcDict.subDict("freezing");
        Tlf_ = freezeDict.lookupOrDefault<scalar>("T_lowerBound", 303.15);
        Tuf_ = freezeDict.lookupOrDefault<scalar>("T_upperBound", 313.15);
        Lf_ = freezeDict.lookupOrDefault<scalar>("latentHeat", 163000.0);
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
    }

    if (pcDict.found("density"))
    {
        const dictionary& densDict = pcDict.subDict("density");
        densityModel_ = densDict.lookupOrDefault<word>("model", "linear");
        rhoRef_ = densDict.lookupOrDefault<scalar>("rhoRef", 1967.0);
        rhoSolid_ = densDict.lookupOrDefault<scalar>("rhoSolid", 1967.0);
        rhoLiquid_ = densDict.lookupOrDefault<scalar>("rhoLiquid", 1850.0);
    }

    if (pcDict.found("thermophysical"))
    {
        const dictionary& thermoDict = pcDict.subDict("thermophysical");
        thermoMode_ = thermoDict.lookupOrDefault<word>("mode", "thermo");

        Cps_ = thermoDict.lookupOrDefault<scalar>("CpSolid", 1980.0);
        Cpl_ = thermoDict.lookupOrDefault<scalar>("CpLiquid", 2320.0);
        ks_ = thermoDict.lookupOrDefault<scalar>("kSolid", 0.50);
        kl_ = thermoDict.lookupOrDefault<scalar>("kLiquid", 0.47);
    }

    if (hysteresisActive_ && mag(Cps_ - Cpl_) > 1e-6)
    {
        Info<< "    EHC PCM Info: Unequal Cp (Cps=" << Cps_ << ", Cpl=" << Cpl_
            << " J/(kg.K)) under stateful hysteresis produces path-dependent sensible enthalpy loops across closed thermal cycles." << endl;
    }

    Info<< "    EHC PCM Parameters loaded for region " << mesh_.name() << ":" << nl
        << "      Melting Range: [" << Tlm_ << " - " << Tum_ << "] K, Latent Heat: " << Lm_ << " J/kg" << nl
        << "      Freezing Range: [" << Tlf_ << " - " << Tuf_ << "] K, Latent Heat: " << Lf_ << " J/kg" << nl
        << "      Hysteresis Active: " << (hysteresisActive_ ? "true" : "false") << nl
        << "      Density Model: " << densityModel_ << " (rhoS=" << rhoSolid_ << ", rhoL=" << rhoLiquid_ << ")" << endl;
}

Foam::scalar Foam::pcmEhcModel::integralAlphaMelt(scalar T, scalar alpha_old) const
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

Foam::scalar Foam::pcmEhcModel::integralAlphaFreeze(scalar T, scalar alpha_old) const
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

void Foam::pcmEhcModel::updateHistory()
{
    liquidFraction_old_ = liquidFraction_;
    heatingTrajectory_old_ = heatingTrajectory_;
}


void Foam::pcmEhcModel::correct()
{
    const volScalarField& T = thermo_.T();
    const volScalarField& Told = thermo_.T().oldTime();
    tmp<volScalarField> tCpThermo = thermo_.Cp();
    tmp<volScalarField> tKappaThermo = thermo_.kappa();
    const volScalarField& CpField = tCpThermo();
    const volScalarField& KappaField = tKappaThermo();

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
            if (deltaTStep > 1e-6)
            {
                isHeating = true;
                heatingTrajectory_[cellI] = 1.0;
            }
            else if (deltaTStep < -1e-6)
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

        scalar alphaL_prev = liquidFraction_old_[cellI];
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

        liquidFraction_[cellI] = alphaL_curr;

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
        if (mag(deltaTStep) > 1e-6)
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
        if (mag(deltaTStep) > 1e-6)
        {
            deltaCp = L * (alphaL_curr - alphaL_prev) / deltaTStep;
        }
        else
        {
            if (isHeating && Tcell > Tlm_ && Tcell < Tum_)
            {
                deltaCp = L / (Tum_ - Tlm_);
            }
            else if (!isHeating && Tcell > Tlf_ && Tcell < Tuf_)
            {
                deltaCp = L / (Tuf_ - Tlf_);
            }
        }

        if (deltaCp < -1e-6)
        {
            WarningInFunction
                << "Negative apparent heat capacity encountered (" << deltaCp
                << " J/(kg.K)) at cell " << cellI << ". Clipping to zero." << endl;
        }

        CpEff_[cellI] = cpBase + max(0.0, deltaCp);

        if (densityModel_ == "constant")
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
        fvPatchScalarField& pLiquidFraction = liquidFraction_.boundaryFieldRef()[patchi];
        const fvPatchScalarField& pLiquidFractionOld = liquidFraction_old_.boundaryField()[patchi];
        fvPatchScalarField& pPhaseState = phaseState_.boundaryFieldRef()[patchi];
        fvPatchScalarField& pCpEff = CpEff_.boundaryFieldRef()[patchi];
        fvPatchScalarField& pRho = rho_.boundaryFieldRef()[patchi];
        fvPatchScalarField& pK = k_.boundaryFieldRef()[patchi];

        const fvPatchScalarField& pCpThermo = CpField.boundaryField()[patchi];
        const fvPatchScalarField& pKappaThermo = KappaField.boundaryField()[patchi];

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
                if (deltaTStep > 1e-6)
                {
                    isHeating = true;
                    pHeatTraj[facei] = 1.0;
                }
                else if (deltaTStep < -1e-6)
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

            scalar alphaL_prev = pLiquidFractionOld[facei];
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

            pLiquidFraction[facei] = alphaL_curr;

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
            if (mag(deltaTStep) > 1e-6)
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
            if (mag(deltaTStep) > 1e-6)
            {
                deltaCp = L * (alphaL_curr - alphaL_prev) / deltaTStep;
            }
            else
            {
                if (isHeating && Tface > Tlm_ && Tface < Tum_)
                {
                    deltaCp = L / (Tum_ - Tlm_);
                }
                else if (!isHeating && Tface > Tlf_ && Tface < Tuf_)
                {
                    deltaCp = L / (Tuf_ - Tlf_);
                }
            }

            if (deltaCp < -1e-6)
            {
                WarningInFunction
                    << "Negative apparent heat capacity encountered (" << deltaCp
                    << " J/(kg.K)) at patch face " << facei << ". Clipping to zero." << endl;
            }

            pCpEff[facei] = cpBase + max(0.0, deltaCp);

            if (densityModel_ == "constant")
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

    liquidFraction_.correctBoundaryConditions();
    phaseState_.correctBoundaryConditions();
    CpEff_.correctBoundaryConditions();
    rho_.correctBoundaryConditions();
    k_.correctBoundaryConditions();
}

// ************************************************************************* //
