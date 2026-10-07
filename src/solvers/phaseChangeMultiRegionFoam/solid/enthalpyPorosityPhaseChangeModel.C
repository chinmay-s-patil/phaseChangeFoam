/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Website:  www.openfoam.com
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

#include "enthalpyPorosityPhaseChangeModel.H"
#include "calculatedFvPatchFields.H"
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
    const basicThermo& thermo
)
:
    phaseChangeModel(mesh, thermo, false),
    lambda_(1.0),
    L_(1.0e5),
    Ts_(300.0),
    Tl_(310.0),
    Cu_(1.0e5),
    q_(0.001),
    Cp_(1000.0),
    kSolid_(1.0),
    kLiquid_(1.0)
{
    readEPDict();

    // Ensure phaseFraction_ has a stored oldTime level
    phaseFraction_.oldTime();

    // Initialize phase fraction from T if not restored from restart
    if (!phaseFractionRestored_)
    {
        const scalarField& T = thermo_.T().primitiveField();
        scalarField& f = phaseFraction_.primitiveFieldRef();
        const scalar window = Tl_ - Ts_;

        forAll(f, i)
        {
            if (window > SMALL)
            {
                f[i] = min(max((T[i] - Ts_)/window, 0.0), 1.0);
            }
            else
            {
                f[i] = (T[i] >= Ts_ ? 1.0 : 0.0);
            }
        }
        phaseFraction_.correctBoundaryConditions();
    }

    // Validate calculated patch types
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
    checkCalculatedPatches(phaseState_);
}



void Foam::enthalpyPorosityPhaseChangeModel::readEPDict()
{
    IOobject dictIO
    (
        "phaseChangeDict",
        mesh_.time().constant(),
        mesh_,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    );

    if (!dictIO.typeHeaderOk<dictionary>(true))
    {
        return;
    }

    IOdictionary phaseChangeDict(dictIO);
    const dictionary& pcDict = phaseChangeDict.subDict("phaseChange");

    // 1. Allowed keys check to catch typos
    for (const word& key : pcDict.toc())
    {
        if
        (
            key != "type"
         && key != "active"
         && key != "forward"
         && key != "porosity"
         && key != "lambda"
         && key != "thermophysical"
         && key != "density"
         && key != "direction"
         && key != "convection"
        )
        {
            FatalIOErrorInFunction(pcDict)
                << "Unknown key '" << key << "' in phaseChange dictionary for model 'enthalpyPorosity'.\n"
                << "Allowed keys are: type, active, forward, porosity, lambda, thermophysical, density, direction, convection."
                << exit(FatalIOError);
        }
    }

    // 2. v1 Restriction: reject hysteresis or direction restrictions
    if (pcDict.found("hysteresis"))
    {
        FatalIOErrorInFunction(pcDict)
            << "Model 'enthalpyPorosity' does not support thermal hysteresis block."
            << exit(FatalIOError);
    }

    if (pcDict.found("reverse"))
    {
        FatalIOErrorInFunction(pcDict)
            << "Model 'enthalpyPorosity' does not support directional 'reverse' block."
            << exit(FatalIOError);
    }

    if (pcDict.found("direction"))
    {
        word dirStr = pcDict.lookupOrDefault<word>("direction", "both");
        if (dirStr != "both")
        {
            FatalIOErrorInFunction(pcDict)
                << "Model 'enthalpyPorosity' only supports 'direction both'."
                << exit(FatalIOError);
        }
    }

    // 3. Read & validate required forward window parameters
    if (!pcDict.found("forward"))
    {
        FatalIOErrorInFunction(pcDict)
            << "Missing required 'forward' subdictionary in phaseChange dictionary for model 'enthalpyPorosity'."
            << exit(FatalIOError);
    }

    const dictionary& fwdDict = pcDict.subDict("forward");
    Ts_ = fwdDict.get<scalar>("T_lowerBound");
    Tl_ = fwdDict.get<scalar>("T_upperBound");
    L_  = fwdDict.get<scalar>("latentHeat");

    if (Tl_ < Ts_)
    {
        FatalIOErrorInFunction(fwdDict)
            << "Invalid temperature window: T_upperBound (liquidus=" << Tl_
            << " K) must be >= T_lowerBound (solidus=" << Ts_ << " K)."
            << exit(FatalIOError);
    }

    if (L_ < 0.0)
    {
        FatalIOErrorInFunction(fwdDict)
            << "Latent heat (latentHeat=" << L_ << " J/kg) cannot be negative."
            << exit(FatalIOError);
    }

    // 4. Read porosity parameters
    if (pcDict.found("porosity"))
    {
        const dictionary& porDict = pcDict.subDict("porosity");
        Cu_ = porDict.lookupOrDefault<scalar>("Cu", porDict.lookupOrDefault<scalar>("A_cu", 1.0e5));
        q_  = porDict.lookupOrDefault<scalar>("q", porDict.lookupOrDefault<scalar>("eps", 0.001));
    }

    // 5. Read & validate relaxation factor lambda
    lambda_ = pcDict.lookupOrDefault<scalar>("lambda", 1.0);
    if (lambda_ <= 0.0 || lambda_ > 1.0)
    {
        FatalIOErrorInFunction(pcDict)
            << "Relaxation factor 'lambda' (" << lambda_ << ") must be in range (0, 1]."
            << exit(FatalIOError);
    }

    // 6. Parallel-safe helper to evaluate background thermo properties across processors
    auto getThermoCp = [&]() -> scalar
    {
        tmp<volScalarField> tCp = thermo_.Cp();
        const scalarField& cpCells = tCp().primitiveField();
        scalar sumCp = 0.0;
        forAll(cpCells, celli)
        {
            sumCp += cpCells[celli];
        }
        label count = cpCells.size();
        reduce(sumCp, sumOp<scalar>());
        reduce(count, sumOp<label>());
        return (count > 0 ? sumCp / count : 1000.0);
    };

    auto getThermoKappa = [&]() -> scalar
    {
        tmp<volScalarField> tK = thermo_.kappa();
        const scalarField& kCells = tK().primitiveField();
        scalar sumK = 0.0;
        forAll(kCells, celli)
        {
            sumK += kCells[celli];
        }
        label count = kCells.size();
        reduce(sumK, sumOp<scalar>());
        reduce(count, sumOp<label>());
        return (count > 0 ? sumK / count : 1.0);
    };

    scalar thermoCpVal = getThermoCp();
    scalar thermoKVal  = getThermoKappa();

    // 7. Read thermophysical parameters & validate single Cp and single rho
    if (pcDict.found("thermophysical"))
    {
        const dictionary& tpDict = pcDict.subDict("thermophysical");
        word mode = tpDict.lookupOrDefault<word>("mode", "thermo");
        if (mode == "custom")
        {
            scalar Cps = tpDict.lookupOrDefault<scalar>("CpSolid", thermoCpVal);
            scalar Cpl = tpDict.lookupOrDefault<scalar>("CpLiquid", thermoCpVal);
            if (mag(Cps - Cpl) > SMALL)
            {
                FatalIOErrorInFunction(tpDict)
                    << "Model 'enthalpyPorosity' v1 requires Cp_solid == Cp_liquid (found CpSolid="
                    << Cps << ", CpLiquid=" << Cpl << ")."
                    << exit(FatalIOError);
            }

            if (mag(Cps - thermoCpVal) / (thermoCpVal + SMALL) > 1e-3)
            {
                FatalIOErrorInFunction(tpDict)
                    << "Custom Cp (" << Cps << " J/(kg K)) disagrees with background thermo Cp ("
                    << thermoCpVal << " J/(kg K)).\n"
                    << "Enthalpy-porosity method requires thermo dict Cp to match custom Cp for energy equation consistency."
                    << exit(FatalIOError);
            }

            Cp_ = Cps;
            kSolid_  = tpDict.lookupOrDefault<scalar>("kSolid", thermoKVal);
            kLiquid_ = tpDict.lookupOrDefault<scalar>("kLiquid", thermoKVal);
        }
        else
        {
            Cp_ = thermoCpVal;
            kSolid_  = thermoKVal;
            kLiquid_ = kSolid_;
        }
    }
    else
    {
        Cp_ = thermoCpVal;
        kSolid_  = thermoKVal;
        kLiquid_ = kSolid_;
    }

    // 8. Validate density
    if (pcDict.found("density"))
    {
        const dictionary& denDict = pcDict.subDict("density");
        if (denDict.found("rhoSolid") && denDict.found("rhoLiquid"))
        {
            scalar rhoS = denDict.lookupOrDefault<scalar>("rhoSolid", 1000.0);
            scalar rhoL = denDict.lookupOrDefault<scalar>("rhoLiquid", 1000.0);
            if (mag(rhoS - rhoL) > SMALL)
            {
                FatalIOErrorInFunction(denDict)
                    << "Model 'enthalpyPorosity' v1 requires rho_solid == rho_liquid (found rhoSolid="
                    << rhoS << ", rhoLiquid=" << rhoL << ")."
                    << exit(FatalIOError);
            }
        }
    }

    // 9. Validate that thermo uses constant Cp (hConst or eConst)
    IOobject tpIO
    (
        "thermophysicalProperties",
        mesh_.time().constant(),
        mesh_,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    );

    if (tpIO.typeHeaderOk<dictionary>(true))
    {
        IOdictionary tpDict(tpIO);
        if (tpDict.found("thermoType"))
        {
            word thermoModel = "";
            if (tpDict.isDict("thermoType"))
            {
                thermoModel = tpDict.subDict("thermoType").lookupOrDefault<word>("thermo", "hConst");
            }
            if (thermoModel != "hConst" && thermoModel != "eConst" && thermoModel != "")
            {
                FatalIOErrorInFunction(pcDict)
                    << "Model 'enthalpyPorosity' requires constant Cp (hConst). Found thermo model '"
                    << thermoModel << "'."
                    << exit(FatalIOError);
            }
        }
    }

    tmp<volScalarField> tCpCheck = thermo_.Cp();
    const scalarField& cpCellsCheck = tCpCheck().primitiveField();
    if (cpCellsCheck.size() > 0)
    {
        scalar minCp = cpCellsCheck[0];
        scalar maxCp = cpCellsCheck[0];
        forAll(cpCellsCheck, i)
        {
            minCp = min(minCp, cpCellsCheck[i]);
            maxCp = max(maxCp, cpCellsCheck[i]);
        }
        reduce(minCp, minOp<scalar>());
        reduce(maxCp, maxOp<scalar>());
        if ((maxCp - minCp) / (minCp + SMALL) > 1e-4)
        {
            FatalIOErrorInFunction(pcDict)
                << "Model 'enthalpyPorosity' requires constant Cp (hConst). Found non-uniform Cp range: ["
                << minCp << ", " << maxCp << "] J/(kg K)."
                << exit(FatalIOError);
        }
    }

    // 10. Validate ratio rhoEff*CpEff / (rho*Cp) == 1 once at construction
    {
        tmp<volScalarField> trhoEff = rhoEff();
        tmp<volScalarField> tCpEff = CpEff();
        tmp<volScalarField> trhoThermo = thermo_.rho();
        tmp<volScalarField> tCpThermo = thermo_.Cp();

        const scalar maxDev = max(mag(trhoEff()*tCpEff()/(trhoThermo()*tCpThermo()) - 1.0)).value();
        if (maxDev > 1e-5)
        {
            FatalIOErrorInFunction(pcDict)
                << "Ratio rhoEff*CpEff/(rho*Cp) = " << maxDev + 1.0
                << " deviates from 1 for model 'enthalpyPorosity'. "
                << "Max deviation: " << maxDev
                << exit(FatalIOError);
        }
    }
}


void Foam::enthalpyPorosityPhaseChangeModel::correct()
{
    if (!active_)
    {
        return;
    }

    // Sync base phaseChangeModel::Cp_ and phaseChangeModel::rho_ fields with thermo
    phaseChangeModel::Cp_.primitiveFieldRef() = thermo_.Cp()().primitiveField();
    phaseChangeModel::Cp_.correctBoundaryConditions();

    phaseChangeModel::rho_.primitiveFieldRef() = thermo_.rho()().primitiveField();
    phaseChangeModel::rho_.correctBoundaryConditions();

    // Update k field
    scalarField& kCells = k_.primitiveFieldRef();
    const scalarField& fCells = phaseFraction_.primitiveField();
    forAll(kCells, i)
    {
        kCells[i] = fCells[i]*kLiquid_ + (1.0 - fCells[i])*kSolid_;
    }
    k_.correctBoundaryConditions();

    #ifdef FULLDEBUG
    // Debug assertion under FULLDEBUG
    tmp<volScalarField> trhoEff = rhoEff();
    tmp<volScalarField> tCpEff = CpEff();
    tmp<volScalarField> trhoThermo = thermo_.rho();
    tmp<volScalarField> tCpThermo = thermo_.Cp();

    const scalar maxDev = max(mag(trhoEff()*tCpEff()/(trhoThermo()*tCpThermo()) - 1.0)).value();
    if (maxDev > 1e-5)
    {
        FatalErrorInFunction
            << "Ratio rhoEff*CpEff/(rho*Cp) = " << maxDev + 1.0
            << " deviates from 1 for model 'enthalpyPorosity'. "
            << "Max deviation: " << maxDev
            << exit(FatalError);
    }
    #endif
}


Foam::tmp<Foam::volScalarField> Foam::enthalpyPorosityPhaseChangeModel::rhoEff() const
{
    return thermo_.rho();
}


Foam::tmp<Foam::volScalarField> Foam::enthalpyPorosityPhaseChangeModel::CpEff() const
{
    return thermo_.Cp();
}


Foam::tmp<Foam::volScalarField> Foam::enthalpyPorosityPhaseChangeModel::latentHeatSource() const
{
    auto tSu = tmp<volScalarField>::New
    (
        IOobject
        (
            "latentHeatSource",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimEnergy/dimVolume/dimTime, 0.0)
    );

    if (!active_ || L_ <= 0.0)
    {
        return tSu;
    }

    tmp<volScalarField> trho = thermo_.rho();
    const volScalarField& rho = trho();
    const surfaceScalarField* phiPtr = mesh_.findObject<surfaceScalarField>("phi");


    dimensionedScalar Ldim("L", dimEnergy/dimMass, L_);

    if (phiPtr)
    {
        tSu.ref() = Ldim * (fvc::ddt(rho, phaseFraction_) + fvc::div(*phiPtr, phaseFraction_, "div(phi,phaseFraction)"));
    }
    else
    {
        tSu.ref() = Ldim * fvc::ddt(rho, phaseFraction_);
    }

    return tSu;
}



Foam::tmp<Foam::volScalarField> Foam::enthalpyPorosityPhaseChangeModel::latentHeatSp() const
{
    return tmp<volScalarField>::New
    (
        IOobject
        (
            "latentHeatSp",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimMass/dimVolume/dimTime, 0.0)
    );
}


Foam::tmp<Foam::volScalarField> Foam::enthalpyPorosityPhaseChangeModel::momentumSp() const
{
    auto tCoeff = tmp<volScalarField>::New
    (
        IOobject
        (
            "momentumSp",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimMass/dimVolume/dimTime, 0.0)
    );

    if (!active_ || Cu_ <= 0.0)
    {
        return tCoeff;
    }

    volScalarField& coeff = tCoeff.ref();
    const scalarField& fCells = phaseFraction_.primitiveField();
    scalarField& coeffCells = coeff.primitiveFieldRef();

    forAll(coeffCells, celli)
    {
        scalar f = fCells[celli];
        coeffCells[celli] = Cu_ * sqr(1.0 - f) / (pow3(f) + q_);
    }

    return tCoeff;
}


Foam::tmp<Foam::volScalarField> Foam::enthalpyPorosityPhaseChangeModel::k() const
{
    return k_;
}


void Foam::enthalpyPorosityPhaseChangeModel::postEnergySolve(volScalarField& he)
{
    if (!active_)
    {
        return;
    }

    const scalarField& T = thermo_.T().primitiveField();
    scalarField& f = phaseFraction_.primitiveFieldRef();
    scalarField& h = he.primitiveFieldRef();

    const scalar window = Tl_ - Ts_;
    const scalar denom = (window > SMALL ? (L_/Cp_ + window) : (L_/Cp_));

    forAll(f, i)
    {
        const scalar Tf = Ts_ + f[i]*window;
        const scalar fNew = min(max(f[i] + lambda_*(T[i] - Tf)/denom, 0.0), 1.0);
        h[i] -= L_*(fNew - f[i]);
        f[i] = fNew;
    }

    // Update boundary conditions of f from patch T (do not apply enthalpy correction on boundaries)
    forAll(phaseFraction_.boundaryFieldRef(), patchi)
    {
        fvPatchScalarField& pf = phaseFraction_.boundaryFieldRef()[patchi];
        const fvPatchScalarField& pT = thermo_.T().boundaryField()[patchi];
        forAll(pf, facei)
        {
            scalar Tb = pT[facei];
            pf[facei] = (window > SMALL ? min(max((Tb - Ts_)/window, 0.0), 1.0) : (Tb >= Ts_ ? 1.0 : 0.0));
        }
    }


    phaseFraction_.correctBoundaryConditions();
    he.correctBoundaryConditions();

    // Update phaseState_ (0 solid, 1 mushy, 2 liquid)
    scalarField& state = phaseState_.primitiveFieldRef();
    forAll(state, i)
    {
        if (f[i] <= 0.0)
        {
            state[i] = 0.0;
        }
        else if (f[i] >= 1.0)
        {
            state[i] = 2.0;
        }
        else
        {
            state[i] = 1.0;
        }
    }
}

// ************************************************************************* //
