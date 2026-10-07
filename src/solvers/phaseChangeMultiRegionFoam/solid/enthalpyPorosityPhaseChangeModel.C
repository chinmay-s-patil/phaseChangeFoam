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
    if (!phaseChangeDict.found("phaseChange"))
    {
        return;
    }

    const dictionary& pcDict = phaseChangeDict.subDict("phaseChange");

    // v1 Restriction: reject hysteresis or direction restrictions
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
        word dirStr = pcDict.get<word>("direction");
        if (dirStr != "both")
        {
            FatalIOErrorInFunction(pcDict)
                << "Model 'enthalpyPorosity' only supports 'direction both'."
                << exit(FatalIOError);
        }
    }

    // Read forward window parameters
    if (pcDict.found("forward"))
    {
        const dictionary& fwdDict = pcDict.subDict("forward");
        Ts_ = fwdDict.get<scalar>("T_lowerBound");
        Tl_ = fwdDict.get<scalar>("T_upperBound");
        L_  = fwdDict.get<scalar>("latentHeat");
    }

    // Read porosity parameters
    if (pcDict.found("porosity"))
    {
        const dictionary& porDict = pcDict.subDict("porosity");
        Cu_ = porDict.getOrDefault<scalar>("Cu", porDict.lookupOrDefault<scalar>("A_cu", 1.0e5));
        q_  = porDict.getOrDefault<scalar>("q", porDict.lookupOrDefault<scalar>("eps", 0.001));
    }

    // Read relaxation factor lambda
    lambda_ = pcDict.getOrDefault<scalar>("lambda", 1.0);

    // Read thermophysical parameters & validate single Cp and single rho
    if (pcDict.found("thermophysical"))
    {
        const dictionary& tpDict = pcDict.subDict("thermophysical");
        word mode = tpDict.getOrDefault<word>("mode", "thermo");
        if (mode == "custom")
        {
            scalar Cps = tpDict.get<scalar>("CpSolid");
            scalar Cpl = tpDict.get<scalar>("CpLiquid");
            if (mag(Cps - Cpl) > SMALL)
            {
                FatalIOErrorInFunction(tpDict)
                    << "Model 'enthalpyPorosity' v1 requires Cp_solid == Cp_liquid (found CpSolid="
                    << Cps << ", CpLiquid=" << Cpl << ")."
                    << exit(FatalIOError);
            }
            Cp_ = Cps;
            kSolid_  = tpDict.get<scalar>("kSolid");
            kLiquid_ = tpDict.get<scalar>("kLiquid");
        }
        else
        {
            Cp_ = thermo_.Cp()()[0];
            kSolid_  = thermo_.kappa()()[0];
            kLiquid_ = kSolid_;
        }
    }
    else
    {
        Cp_ = thermo_.Cp()()[0];
        kSolid_  = thermo_.kappa()()[0];
        kLiquid_ = kSolid_;
    }

    // Validate density
    if (pcDict.found("density"))
    {
        const dictionary& denDict = pcDict.subDict("density");
        if (denDict.found("rhoSolid") && denDict.found("rhoLiquid"))
        {
            scalar rhoS = denDict.get<scalar>("rhoSolid");
            scalar rhoL = denDict.get<scalar>("rhoLiquid");
            if (mag(rhoS - rhoL) > SMALL)
            {
                FatalIOErrorInFunction(denDict)
                    << "Model 'enthalpyPorosity' v1 requires rho_solid == rho_liquid (found rhoSolid="
                    << rhoS << ", rhoLiquid=" << rhoL << ")."
                    << exit(FatalIOError);
            }
        }
    }
}


void Foam::enthalpyPorosityPhaseChangeModel::correct()
{
    if (!active_)
    {
        return;
    }

    // Update k field
    scalarField& kCells = k_.primitiveFieldRef();
    const scalarField& fCells = phaseFraction_.primitiveField();
    forAll(kCells, i)
    {
        kCells[i] = fCells[i]*kLiquid_ + (1.0 - fCells[i])*kSolid_;
    }
    k_.correctBoundaryConditions();
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

    const auto& rho = mesh_.lookupObject<volScalarField>("rho");
    const auto& phi = mesh_.lookupObject<surfaceScalarField>("phi");

    tSu.ref() = L_ * (fvc::ddt(rho, phaseFraction_) + fvc::div(phi, phaseFraction_, "div(phi,phaseFraction)"));
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
