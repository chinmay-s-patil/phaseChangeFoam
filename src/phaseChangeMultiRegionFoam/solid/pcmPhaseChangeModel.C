/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Website:  www.openfoam.com
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

#include "pcmPhaseChangeModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(pcmPhaseChangeModel, 0);
    defineRunTimeSelectionTable(pcmPhaseChangeModel, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pcmPhaseChangeModel::pcmPhaseChangeModel
(
    const fvMesh& mesh,
    const solidThermo& thermo
)
:
    mesh_(mesh),
    thermo_(thermo),
    active_(false),
    liquidFraction_
    (
        IOobject
        (
            "liquidFraction",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("liquidFraction", dimless, 0.0),
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
    )
{}


// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::pcmPhaseChangeModel> Foam::pcmPhaseChangeModel::New
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

        bool active = phaseChangeDict.lookupOrDefault<bool>("active", false);

        if (active && phaseChangeDict.found("phaseChange"))
        {
            const dictionary& pcDict = phaseChangeDict.subDict("phaseChange");
            bool pcActive = pcDict.lookupOrDefault<bool>("active", false);

            if (pcActive)
            {
                word modeName = pcDict.lookupOrDefault<word>("phaseChangeMode", "EHC");

                Info<< "    Phase change ACTIVE for region " << mesh.name()
                    << " using mode: " << modeName << endl;

                auto cstrIter = dictionaryConstructorTablePtr_->cfind(modeName);

                if (!cstrIter.good())
                {
                    FatalErrorInFunction
                        << "Unknown phaseChangeMode " << modeName
                        << " for region " << mesh.name() << nl << nl
                        << "Valid options are: "
                        << dictionaryConstructorTablePtr_->toc()
                        << exit(FatalError);
                }

                return autoPtr<pcmPhaseChangeModel>(cstrIter()(mesh, thermo));
            }
        }
    }

    Info<< "    Phase change INACTIVE for region " << mesh.name() << endl;

    // Return dummy inactive model if not configured
    auto cstrIter = dictionaryConstructorTablePtr_->cfind("none");

    if (cstrIter.good())
    {
        return autoPtr<pcmPhaseChangeModel>(cstrIter()(mesh, thermo));
    }

    // Default fallback: return nullptr if inactive
    return autoPtr<pcmPhaseChangeModel>(nullptr);
}

// ************************************************************************* //
