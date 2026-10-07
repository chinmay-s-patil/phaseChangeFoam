    CombustionModel<rhoReactionThermo>& reaction = reactionFluid[i];

    rhoReactionThermo& thermo = reaction.thermo();
    thermo.validate(args.executable(), "h", "e");

    basicMultiComponentMixture& composition = thermo.composition();
    PtrList<volScalarField>& Y = composition.Y();

    label inertIndex = -1;
    if (Y.size())
    {
        const word inertSpecie(thermo.get<word>("inertSpecie"));
        inertIndex = composition.species().find(inertSpecie);
        if (inertIndex < 0)
        {
            FatalIOErrorIn(args.executable().c_str(), thermo)
                << "Inert specie " << inertSpecie
                << " not found in available species "
                << flatOutput(composition.species())
                << exit(FatalIOError);
        }
    }

    volScalarField& rho = rhoFluid[i];
    volVectorField& U = UFluid[i];
    surfaceScalarField& phi = phiFluid[i];

    compressible::turbulenceModel& turb = turbulence[i];

    volScalarField& p = thermo.p();
    const volScalarField& psi = thermo.psi();

    volScalarField& p_rgh = p_rghFluid[i];

    const volScalarField& gh = ghFluid[i];
    const surfaceScalarField& ghf = ghfFluid[i];

    multivariateSurfaceInterpolationScheme<scalar>::fieldTable& fields =
        fieldsFluid[i];

    volScalarField& Qdot = QdotFluid[i];

    radiation::radiationModel& rad = radiation[i];

    fluidPhaseChangeModel& phaseChange = fluidPhaseChangeModels[i];
    phaseChangeModel& pcmModel = pcmModelsFluid[i];

    IOMRFZoneList& MRF = MRFfluid[i];
    fv::options& fvOptions = fluidFvOptions[i];

    fvVectorMatrix& UEqn = UEqFluid[i];

    const dimensionedScalar initialMass
    (
        "initialMass",
        dimMass,
        initialMassFluid[i]
    );

    bool frozenFlow = frozenFlowFluid[i];

    const label pRefCell = pRefCellFluid[i];
    const scalar pRefValue = pRefValueFluid[i];

