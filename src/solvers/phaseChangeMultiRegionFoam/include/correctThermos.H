forAll(fluidRegions, i)
{
    rhoThermo& thermo = thermoFluid[i];
    radiation::radiationModel& rad = radiation[i];
    fv::options& fvOptions = fluidFvOptions[i];
    volScalarField& he = thermo.he();
    fvOptions.correct(he);
    thermo.correct();
    phaseChangeModel& pcmModel = pcmModelsFluid[i];
    if (pcmModel.active())
    {
        pcmModel.postEnergySolve(he);
        thermo.correct();
    }
    rad.correct();
}

forAll(solidRegions, i)
{
    solidThermo& thermo = thermos[i];
    fv::options& fvOptions = solidHeatSources[i];
    volScalarField& h = thermo.he();
    fvOptions.correct(h);
    thermo.correct();
    phaseChangeModel& pcModel = phaseChangeModels[i];
    if (pcModel.active())
    {
        pcModel.postEnergySolve(h);
        thermo.correct();
    }
}


