/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Website:  www.openfoam.com
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

#include "phaseChangeModel.H"
#include <cctype>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(phaseChangeModel, 0);
    defineRunTimeSelectionTable(phaseChangeModel, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseChangeModel::phaseChangeModel
(
    const fvMesh& mesh,
    const solidThermo& thermo
)
:
    mesh_(mesh),
    thermo_(thermo),
    active_(false),
    phaseFraction_
    (
        IOobject
        (
            "phaseFraction",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("phaseFraction", dimless, 0.0),
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
    ),
    phaseFractionRestored_(false),
    suppressConvection_(true)
{
    phaseFractionRestored_ = phaseFraction_.headerOk();

    // Backward compatibility for disk reads: if liquidFraction exists on disk but phaseFraction doesn't
    IOobject liquidFractionIO
    (
        "liquidFraction",
        mesh.time().timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    );

    if (!phaseFraction_.headerOk() && liquidFractionIO.typeHeaderOk<volScalarField>(true))
    {
        volScalarField lFrac(liquidFractionIO, mesh);
        phaseFraction_ == lFrac;
        phaseFractionRestored_ = true;
    }
}


void Foam::phaseChangeModel::readConvectionDict(const dictionary& pcDict)
{
    if (pcDict.found("convection"))
    {
        suppressConvection_ =
            pcDict.subDict("convection").getOrDefault<bool>("suppress", true);
    }

    if (!suppressConvection_)
    {
        FatalIOErrorInFunction(pcDict)
            << "phaseChange.convection.suppress = false requested for region "
            << mesh_.name() << ", but phase change is solved in solid-type "
            << "regions without a momentum equation; convection cannot be "
            << "resolved. Use suppress true (conduction-only)."
            << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::phaseChangeModel> Foam::phaseChangeModel::New
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

        bool active = phaseChangeDict.lookupOrDefault<bool>("active", true);

        if (active && phaseChangeDict.found("phaseChange"))
        {
            const dictionary& pcDict = phaseChangeDict.subDict("phaseChange");
            bool pcActive = pcDict.lookupOrDefault<bool>("active", true);

            if (pcActive)
            {
                word modeName = pcDict.lookupOrDefault<word>("type", pcDict.lookupOrDefault<word>("phaseChangeMode", "ehc"));

                Info<< "    Phase change ACTIVE for region " << mesh.name()
                    << " using mode: " << modeName << endl;

                auto cstrIter = dictionaryConstructorTablePtr_->cfind(modeName);

                if (!cstrIter.good())
                {
                    auto toLowerStr = [](const word& w)
                    {
                        std::string s(w);
                        for (char& c : s)
                        {
                            c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
                        }
                        return s;
                    };

                    for (auto iter = dictionaryConstructorTablePtr_->cbegin(); iter != dictionaryConstructorTablePtr_->cend(); ++iter)
                    {
                        if (toLowerStr(iter.key()) == toLowerStr(modeName))
                        {
                            cstrIter = iter;
                            break;
                        }
                    }
                }

                if (!cstrIter.good())
                {
                    FatalErrorInFunction
                        << "Unknown phaseChange mode " << modeName
                        << " for region " << mesh.name() << nl << nl
                        << "Valid options are: "
                        << dictionaryConstructorTablePtr_->toc()
                        << exit(FatalError);
                }

                return autoPtr<phaseChangeModel>(cstrIter()(mesh, thermo));
            }
        }
    }

    Info<< "    Phase change INACTIVE for region " << mesh.name() << endl;

    // Return dummy inactive model if not configured
    auto cstrIter = dictionaryConstructorTablePtr_->cfind("none");

    if (cstrIter.good())
    {
        return autoPtr<phaseChangeModel>(cstrIter()(mesh, thermo));
    }

    // Default fallback: return nullptr if inactive
    return autoPtr<phaseChangeModel>(nullptr);
}

// ************************************************************************* //
