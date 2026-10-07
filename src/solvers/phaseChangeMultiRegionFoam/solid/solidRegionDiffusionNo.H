scalar DiNum = -GREAT;

forAll(solidRegions, i)
{
    //- Note: do not use setRegionSolidFields.H to avoid double registering Cp
    //#include "setRegionSolidFields.H"
    const solidThermo& thermo = thermos[i];
    const phaseChangeModel& pcModel = phaseChangeModels[i];

    tmp<volScalarField> magKappa =
        pcModel.active()
      ? pcModel.k()
      : (thermo.isotropic() ? thermo.kappa() : tmp<volScalarField>(mag(thermo.Kappa())));

    tmp<volScalarField> tcp = pcModel.active() ? pcModel.CpEff() : tmp<volScalarField>(thermo.Cp());
    tmp<volScalarField> trho = pcModel.active() ? pcModel.rhoEff() : tmp<volScalarField>(thermo.rho());

    volScalarField Cprho("Cprho", trho() * tcp());

    DiNum = max
    (
        solidRegionDiffNo
        (
            solidRegions[i],
            runTime,
            Cprho,
            magKappa()
        ),
        DiNum
    );

}
