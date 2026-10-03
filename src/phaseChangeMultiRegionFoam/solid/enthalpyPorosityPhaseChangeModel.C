/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Website:  www.openfoam.com
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

#include "enthalpyPorosityPhaseChangeModel.H"
#include "addToRunTimeSelectionTable.H"
#include "calculatedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(enthalpyPorosityPhaseChangeModel, 0);
    addToRunTimeSelectionTable(phaseChangeModel, enthalpyPorosityPhaseChangeModel, dictionary);
}

// * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

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
    Tlf_(303.15),
    Tuf_(313.15),
    Lf_(163000.0),
    hysteresisActive_(true),
    reversalTol_(1e-6),
    beta_(50e-6),
    Cu_(1e5),
    q_(1e-2),
    suppressConvection_(false),
    densityModel_("thermo"),
    rhoRef_(0.0),
    rhoSolid_(0.0),
    rhoLiquid_(0.0),
    thermoMode_("thermo"),
    Cps_(1980.0),
    Cpl_(2320.0),
    ks_(0.50),
    kl_(0.47),
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

    // Validate that all non-constraint patches are 'calculated'
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
    checkCalculatedPatches(Cp_);
    checkCalculatedPatches(rho_);
    checkCalculatedPatches(k_);
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
    suppressConvection_ = suppressConvection();

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

    if (pcDict.found("hysteresis"))
    {
        const dictionary& hysDict = pcDict.subDict("hysteresis");
        hysteresisActive_ = hysDict.lookupOrDefault<bool>("active", true);
        reversalTol_ = hysDict.lookupOrDefault<scalar>("reversalTol", 1e-6);
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
}

void Foam::enthalpyPorosityPhaseChangeModel::correct()
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

        scalar Tl = isHeating ? Tlm_ : Tlf_;
        scalar Tu = isHeating ? Tum_ : Tuf_;

        scalar alpha_m_curr = (Tcell <= Tlm_) ? 0.0 : ((Tcell >= Tum_) ? 1.0 : (Tcell - Tlm_) / (Tum_ - Tlm_));
        scalar alpha_f_curr = (Tcell <= Tlf_) ? 0.0 : ((Tcell >= Tuf_) ? 1.0 : (Tcell - Tlf_) / (Tuf_ - Tlf_));

        scalar alphaL_prev = phaseFraction_old_[cellI];
        if (!hysteresisActive_)
        {
            alphaL_prev = (ToldCell <= Tl) ? 0.0 : ((ToldCell >= Tu) ? 1.0 : (ToldCell - Tl) / (Tu - Tl));
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

        Cp_[cellI] = (1.0 - alphaL_curr) * cpsVal + alphaL_curr * cplVal;

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

    forAll(T.boundaryField(), patchi)
    {
        const fvPatchScalarField& pT = T.boundaryField()[patchi];
        const fvPatchScalarField& pTold = Told.boundaryField()[patchi];
        fvPatchScalarField& pHeatTraj = heatingTrajectory_.boundaryFieldRef()[patchi];
        const fvPatchScalarField& pHeatTrajOld = heatingTrajectory_old_.boundaryField()[patchi];
        fvPatchScalarField& pPhaseFraction = phaseFraction_.boundaryFieldRef()[patchi];
        const fvPatchScalarField& pPhaseFractionOld = phaseFraction_old_.boundaryField()[patchi];
        fvPatchScalarField& pPhaseState = phaseState_.boundaryFieldRef()[patchi];
        fvPatchScalarField& pCp = Cp_.boundaryFieldRef()[patchi];
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

            scalar Tl = isHeating ? Tlm_ : Tlf_;
            scalar Tu = isHeating ? Tum_ : Tuf_;

            scalar alpha_m_curr = (Tface <= Tlm_) ? 0.0 : ((Tface >= Tum_) ? 1.0 : (Tface - Tlm_) / (Tum_ - Tlm_));
            scalar alpha_f_curr = (Tface <= Tlf_) ? 0.0 : ((Tface >= Tuf_) ? 1.0 : (Tface - Tlf_) / (Tuf_ - Tlf_));

            scalar alphaL_prev = pPhaseFractionOld[facei];
            if (!hysteresisActive_)
            {
                alphaL_prev = (ToldFace <= Tl) ? 0.0 : ((ToldFace >= Tu) ? 1.0 : (ToldFace - Tl) / (Tu - Tl));
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

            pCp[facei] = (1.0 - alphaL_curr) * cpsVal + alphaL_curr * cplVal;

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

    Cp_.correctBoundaryConditions();
    rho_.correctBoundaryConditions();
    k_.correctBoundaryConditions();
}

void Foam::enthalpyPorosityPhaseChangeModel::updateHistory()
{
    phaseFraction_old_ = phaseFraction_;
    heatingTrajectory_old_ = heatingTrajectory_;
}

Foam::tmp<Foam::volScalarField> Foam::enthalpyPorosityPhaseChangeModel::latentHeatSource() const
{
    tmp<volScalarField> tSu
    (
        volScalarField::New
        (
            "latentHeatSource",
            mesh_,
            dimensionedScalar("zero", dimEnergy/dimTime/dimVolume, 0.0)
        )
    );

    scalarField& Su = tSu.ref().primitiveFieldRef();
    const scalarField& T = thermo_.T().primitiveField();
    const scalarField& Told = thermo_.T().oldTime().primitiveField();
    const scalarField& h = thermo_.he().primitiveField();
    const scalar rDt = 1.0/mesh_.time().deltaTValue();

    forAll(T, i)
    {
        bool isHeating = (heatingTrajectory_[i] > 0.5);
        scalar Tl = isHeating ? Tlm_ : Tlf_;
        scalar Tu = isHeating ? Tum_ : Tuf_;
        scalar L = isHeating ? Lm_ : Lf_;
        scalar dT_mush = Tu - Tl;

        const scalar rhoL = rho_[i]*L*rDt;
        const scalar aCurr = phaseFraction_[i];
        const scalar aOld = phaseFraction_old_[i];

        if (T[i] <= Tl && aOld <= 0.0)
        {
            Su[i] = 0.0;
        }
        else if (T[i] >= Tu && aOld >= 1.0)
        {
            Su[i] = 0.0;
        }
        else
        {
            scalar dAlpha = 1.0 / dT_mush;
            if (T[i] <= Tl || T[i] >= Tu)
            {
                scalar deltaT = mag(T[i] - Told[i]);
                if (deltaT > dT_mush && mag(aCurr - aOld) > 1e-6)
                {
                    dAlpha = mag(aCurr - aOld) / deltaT;
                }
            }
            const scalar Sp = rhoL * dAlpha / max(Cp_[i], 1e-10);
            Su[i] = rhoL * (aCurr - aOld) - Sp * h[i];
        }
    }
    return tSu;
}

Foam::tmp<Foam::volScalarField> Foam::enthalpyPorosityPhaseChangeModel::latentHeatSp() const
{
    tmp<volScalarField> tSp
    (
        volScalarField::New
        (
            "latentHeatSp",
            mesh_,
            dimensionedScalar("zero", dimMass/dimVolume/dimTime, 0.0)
        )
    );

    scalarField& Sp = tSp.ref().primitiveFieldRef();
    const scalarField& T = thermo_.T().primitiveField();
    const scalarField& Told = thermo_.T().oldTime().primitiveField();
    const scalar rDt = 1.0/mesh_.time().deltaTValue();

    forAll(T, i)
    {
        bool isHeating = (heatingTrajectory_[i] > 0.5);
        scalar Tl = isHeating ? Tlm_ : Tlf_;
        scalar Tu = isHeating ? Tum_ : Tuf_;
        scalar L = isHeating ? Lm_ : Lf_;
        scalar dT_mush = Tu - Tl;

        const scalar aCurr = phaseFraction_[i];
        const scalar aOld = phaseFraction_old_[i];

        if (T[i] <= Tl && aOld <= 0.0)
        {
            Sp[i] = 0.0;
        }
        else if (T[i] >= Tu && aOld >= 1.0)
        {
            Sp[i] = 0.0;
        }
        else
        {
            scalar dAlpha = 1.0 / dT_mush;
            if (T[i] <= Tl || T[i] >= Tu)
            {
                scalar deltaT = mag(T[i] - Told[i]);
                if (deltaT > dT_mush && mag(aCurr - aOld) > 1e-6)
                {
                    dAlpha = mag(aCurr - aOld) / deltaT;
                }
            }
            Sp[i] = rho_[i] * L * rDt * dAlpha / max(Cp_[i], 1e-10);
        }
    }
    return tSp;
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

    if (!suppressConvection_ && mesh_.foundObject<volVectorField>("U"))
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
