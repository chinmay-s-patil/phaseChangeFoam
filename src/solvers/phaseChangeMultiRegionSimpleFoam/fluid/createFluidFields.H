#include "fluidPhaseChangeModel.H"

// Initialise fluid field pointer lists
PtrList<rhoReactionThermo> thermoFluid(fluidRegions.size());
PtrList<volScalarField> rhoFluid(fluidRegions.size());
PtrList<volVectorField> UFluid(fluidRegions.size());
PtrList<surfaceScalarField> phiFluid(fluidRegions.size());
PtrList<fluidPhaseChangeModel> fluidPhaseChangeModels(fluidRegions.size());
PtrList<phaseChangeModel> pcmModelsFluid(fluidRegions.size());
PtrList<uniformDimensionedScalarField> hRefFluid(fluidRegions.size());
PtrList<volScalarField> ghFluid(fluidRegions.size());
PtrList<surfaceScalarField> ghfFluid(fluidRegions.size());
PtrList<compressible::turbulenceModel> turbulence(fluidRegions.size());
PtrList<CombustionModel<rhoReactionThermo>> reactionFluid(fluidRegions.size());
PtrList<volScalarField> p_rghFluid(fluidRegions.size());
PtrList<radiation::radiationModel> radiation(fluidRegions.size());
PtrList<multivariateSurfaceInterpolationScheme<scalar>::fieldTable>
    fieldsFluid(fluidRegions.size());
PtrList<volScalarField> QdotFluid(fluidRegions.size());

List<scalar> initialMassFluid(fluidRegions.size());
List<label> pRefCellFluid(fluidRegions.size(), -1);
List<scalar> pRefValueFluid(fluidRegions.size(), Zero);
List<bool> frozenFlowFluid(fluidRegions.size(), false);

PtrList<dimensionedScalar> rhoMax(fluidRegions.size());
PtrList<dimensionedScalar> rhoMin(fluidRegions.size());

PtrList<IOMRFZoneList> MRFfluid(fluidRegions.size());
PtrList<fv::options> fluidFvOptions(fluidRegions.size());

PtrList<fvVectorMatrix> UEqFluid(fluidRegions.size());

const uniformDimensionedVectorField& g = meshObjects::gravity::New(runTime);

// Populate fluid field pointer lists
forAll(fluidRegions, i)
{
    Info<< "*** Reading fluid mesh thermophysical properties for region "
        << fluidRegions[i].name() << nl << endl;

    Info<< "    Adding to thermoFluid\n" << endl;

    thermoFluid.set
    (
        i,
        rhoReactionThermo::New(fluidRegions[i]).ptr()
    );

    Info<< "    Adding to rhoFluid\n" << endl;
    rhoFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "rho",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            thermoFluid[i].rho()
        )
    );

    Info<< "    Adding to UFluid\n" << endl;
    UFluid.set
    (
        i,
        new volVectorField
        (
            IOobject
            (
                "U",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i]
        )
    );

    Info<< "    Adding to phiFluid\n" << endl;
    phiFluid.set
    (
        i,
        new surfaceScalarField
        (
            IOobject
            (
                "phi",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            linearInterpolate(rhoFluid[i]*UFluid[i])
                & fluidRegions[i].Sf()
        )
    );

    Info<< "    Adding to hRefFluid\n" << endl;
    hRefFluid.set
    (
        i,
        new uniformDimensionedScalarField
        (
            IOobject
            (
                "hRef",
                runTime.constant(),
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            dimensionedScalar(word::null, dimLength, Zero)
        )
    );

    dimensionedScalar ghRef
    (
        mag(g.value()) > SMALL
      ? g & (cmptMag(g.value())/mag(g.value()))*hRefFluid[i]
      : dimensionedScalar("ghRef", g.dimensions()*dimLength, 0)
    );

    Info<< "    Adding to ghFluid\n" << endl;
    ghFluid.set
    (
        i,
        new volScalarField
        (
            "gh",
            (g & fluidRegions[i].C()) - ghRef
        )
    );

    Info<< "    Adding to ghfFluid\n" << endl;
    ghfFluid.set
    (
        i,
        new surfaceScalarField
        (
            "ghf",
            (g & fluidRegions[i].Cf()) - ghRef
        )
    );

    Info<< "    Adding to turbulence\n" << endl;
    turbulence.set
    (
        i,
        compressible::turbulenceModel::New
        (
            rhoFluid[i],
            UFluid[i],
            phiFluid[i],
            thermoFluid[i]
        ).ptr()
    );

    Info<< "    Adding to reactionFluid\n" << endl;
    reactionFluid.set
    (
        i,
        CombustionModel<rhoReactionThermo>::New
        (
            thermoFluid[i],
            turbulence[i]
        )
    );

    p_rghFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "p_rgh",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i]
        )
    );

    // Force p_rgh to be consistent with p
    p_rghFluid[i] = thermoFluid[i].p() - rhoFluid[i]*ghFluid[i];

    fluidRegions[i].setFluxRequired(p_rghFluid[i].name());

    radiation.set
    (
        i,
        radiation::radiationModel::New(thermoFluid[i].T())
    );

    initialMassFluid[i] = fvc::domainIntegrate(rhoFluid[i]).value();

    Info<< "    Adding to fieldsFluid\n" << endl;
    fieldsFluid.set
    (
        i,
        new multivariateSurfaceInterpolationScheme<scalar>::fieldTable
    );
    forAll(thermoFluid[i].composition().Y(), j)
    {
        fieldsFluid[i].add(thermoFluid[i].composition().Y()[j]);
    }
    fieldsFluid[i].add(thermoFluid[i].he());

    Info<< "    Adding to QdotFluid\n" << endl;
    QdotFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "Qdot",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i],
            dimensionedScalar(dimEnergy/dimVolume/dimTime, Zero)
        )
    );

    const dictionary& simpleDict =
        fluidRegions[i].solutionDict().subDict("SIMPLE");

    setRefCell
    (
        thermoFluid[i].p(),
        p_rghFluid[i],
        simpleDict,
        pRefCellFluid[i],
        pRefValueFluid[i]
    );

    simpleDict.readIfPresent("frozenFlow", frozenFlowFluid[i]);

    rhoMax.set
    (
        i,
        new dimensionedScalar("rhoMax", dimDensity, GREAT, simpleDict)
    );

    rhoMin.set
    (
        i,
        new dimensionedScalar("rhoMin", dimDensity, Zero, simpleDict)
    );

    Info<< "    Adding MRF\n" << endl;
    MRFfluid.set
    (
        i,
        new IOMRFZoneList(fluidRegions[i])
    );

    Info<< "    Adding fvOptions\n" << endl;
    fluidFvOptions.set
    (
        i,
        new fv::options(fluidRegions[i])
    );

    Info<< "    Adding fluid phase change model\n" << endl;
    fluidPhaseChangeModels.set
    (
        i,
        fluidPhaseChangeModel::New
        (
            fluidRegions[i],
            thermoFluid[i],
            UFluid[i],
            phiFluid[i]
        )
    );

    if (fluidPhaseChangeModels[i].active() && Foam::isA<Foam::leeFluidPhaseChangeModel>(fluidPhaseChangeModels[i]))
    {
        const Foam::leeFluidPhaseChangeModel& leeModel =
            Foam::refCast<const Foam::leeFluidPhaseChangeModel>(fluidPhaseChangeModels[i]);

        const Foam::basicMultiComponentMixture& composition = thermoFluid[i].composition();
        if (composition.Y().size())
        {
            const Foam::word inertSpecie(thermoFluid[i].get<Foam::word>("inertSpecie"));
            const Foam::label inertIndex = composition.species().find(inertSpecie);

            if (inertIndex >= 0 && (inertIndex == leeModel.liquidIndex() || inertIndex == leeModel.vaporIndex()))
            {
                FatalErrorInFunction
                    << "inertSpecie '" << inertSpecie
                    << "' in fluid region " << fluidRegions[i].name()
                    << " must not be a phase-changing specie ("
                    << leeModel.liquidName() << " or " << leeModel.vaporName() << ")."
                    << Foam::exit(Foam::FatalError);
            }
        }
    }

    Info<< "    Adding PCM phase change model\n" << endl;
    pcmModelsFluid.set
    (
        i,
        phaseChangeModel::New
        (
            fluidRegions[i],
            thermoFluid[i]
        )
    );

    if (pcmModelsFluid[i].active())
    {
        const phaseChangeModel& pcm = pcmModelsFluid[i];
        const scalar Cps = pcm.Cps();
        const scalar Cpl = pcm.Cpl();
        const scalar rhos = pcm.rhoSolid();
        const scalar rhol = pcm.rhoLiquid();

        tmp<volScalarField> tCpThermo = thermoFluid[i].Cp();
        tmp<volScalarField> trhoThermo = thermoFluid[i].rho();
        const scalarField& cpCells = tCpThermo().primitiveField();
        const scalarField& rhoCells = trhoThermo().primitiveField();

        const scalar maxCpDevS = gMax(mag(cpCells - Cps))/(Cps + 1e-12);
        const scalar maxCpDevL = gMax(mag(cpCells - Cpl))/(Cpl + 1e-12);
        const scalar maxCpDevSL = mag(Cps - Cpl)/(Cps + 1e-12);

        const scalar maxRhoDevS = gMax(mag(rhoCells - rhos))/(rhos + 1e-12);
        const scalar maxRhoDevL = gMax(mag(rhoCells - rhol))/(rhol + 1e-12);
        const scalar maxRhoDevSL = mag(rhos - rhol)/(rhos + 1e-12);

        scalar maxDev = max(maxCpDevS, maxCpDevL);
        maxDev = max(maxDev, maxCpDevSL);
        maxDev = max(maxDev, maxRhoDevS);
        maxDev = max(maxDev, maxRhoDevL);
        maxDev = max(maxDev, maxRhoDevSL);

        if (maxDev > 1e-5)
        {
            FatalErrorInFunction
                << "PCM model '" << pcm.type()
                << "' in fluid region " << fluidRegions[i].name()
                << " requires equal solid/liquid heat capacities and densities matching background thermophysical properties "
                << "(Cp_solid == Cp_liquid == thermo Cp, and rho_solid == rho_liquid == thermo rho). "
                << "Found max deviation = " << maxDev << ".\n"
                << "Parameters: Cp_solid=" << Cps << ", Cp_liquid=" << Cpl
                << ", rho_solid=" << rhos << ", rho_liquid=" << rhol << "."
                << exit(FatalError);
        }
    }

    UEqFluid.set
    (
        i,
        new fvVectorMatrix(UFluid[i], dimForce)
    );

    turbulence[i].validate();
}

