/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
\*---------------------------------------------------------------------------*/

#include "fluidPhaseChangeModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fluidPhaseChangeModel, 0);
    defineRunTimeSelectionTable(fluidPhaseChangeModel, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fluidPhaseChangeModel::fluidPhaseChangeModel
(
    const fvMesh& mesh,
    const rhoReactionThermo& thermo,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    IOdictionary
    (
        IOobject
        (
            "phaseChangeDict",
            mesh.time().constant(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    ),
    mesh_(mesh),
    thermo_(thermo),
    U_(U),
    phi_(phi),
    dict_(subOrEmptyDict("phaseChange")),
    active_(dict_.getOrDefault<bool>("active", false)),
    zeroSpeciesSource_
    (
        IOobject
        (
            "zeroSpeciesSource",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimMass/dimVolume/dimTime, Zero)
    )
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::fluidPhaseChangeModel>
Foam::fluidPhaseChangeModel::New
(
    const fvMesh& mesh,
    const rhoReactionThermo& thermo,
    const volVectorField& U,
    const surfaceScalarField& phi
)
{
    IOobject dictHeader
    (
        "phaseChangeDict",
        mesh.time().constant(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    );

    word modelType("none");

    if (dictHeader.typeHeaderOk<IOdictionary>(true))
    {
        IOdictionary dict(dictHeader);
        if (dict.found("phaseChange"))
        {
            const dictionary& subDict = dict.subDict("phaseChange");
            const bool active = subDict.getOrDefault<bool>("active", false);
            if (active)
            {
                bool typeFound = false;
                if (subDict.found("speciesModel"))
                {
                    modelType = subDict.get<word>("speciesModel");
                    typeFound = true;
                }
                else if (subDict.found("type"))
                {
                    modelType = subDict.get<word>("type");
                    typeFound = true;
                }
                else if (subDict.found("phaseChangeMode"))
                {
                    modelType = subDict.get<word>("phaseChangeMode");
                    typeFound = true;
                }

                auto toLowerStr = [](const word& w)
                {
                    std::string s(w);
                    for (char& c : s)
                    {
                        c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
                    }
                    return s;
                };

                std::string modeLower = toLowerStr(modelType);
                if (modeLower == "ehc" || modeLower == "ehcdarcy" || modeLower == "enthalpyporosity")
                {
                    modelType = "none";
                }

                bool speciesKeysFound = (
                    subDict.found("liquid") ||
                    subDict.found("vapor") ||
                    subDict.found("C_evap") ||
                    subDict.found("C_cond") ||
                    subDict.found("T_sat") ||
                    subDict.found("rate")
                );

                if (!typeFound)
                {
                    if (speciesKeysFound)
                    {
                        FatalIOErrorInFunction(subDict)
                            << "Active phaseChange dictionary in region " << mesh.name()
                            << " contains species phase change parameters (liquid, vapor, C_evap, etc.), "
                            << "but missing 'type' or 'speciesModel' keyword."
                            << exit(FatalIOError);
                    }
                    else if (!subDict.found("meltingModel") && !subDict.found("porosity") && !subDict.found("forward"))
                    {
                        FatalIOErrorInFunction(subDict)
                            << "Active phaseChange dictionary in region " << mesh.name()
                            << " is missing 'type', 'speciesModel', or 'meltingModel' keyword."
                            << exit(FatalIOError);
                    }
                }
                else if (modelType == "none" || modeLower == "none" || modeLower == "off")
                {
                    if (speciesKeysFound)
                    {
                        FatalIOErrorInFunction(subDict)
                            << "Fluid species phase change model resolved to 'none' for region " << mesh.name()
                            << ", but species phase change parameters (liquid, vapor, C_evap, etc.) were found in phaseChangeDict."
                            << exit(FatalIOError);
                    }
                    modelType = "none";
                }
                else
                {
                    if (!dictionaryConstructorTablePtr_->cfind(modelType).good())
                    {
                        for (auto iter = dictionaryConstructorTablePtr_->cbegin(); iter != dictionaryConstructorTablePtr_->cend(); ++iter)
                        {
                            if (toLowerStr(iter.key()) == modeLower)
                            {
                                modelType = iter.key();
                                break;
                            }
                        }
                    }
                }
            }
        }
    }

    Info<< "Selecting fluid phase change model type " << modelType
        << " for region " << mesh.name() << endl;

    auto cstrIter = dictionaryConstructorTablePtr_->cfind(modelType);

    if (!cstrIter.good())
    {
        FatalErrorInFunction
            << "Unknown fluidPhaseChangeModel type "
            << modelType << nl << nl
            << "Valid fluidPhaseChangeModel types are:" << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<fluidPhaseChangeModel>(cstrIter()(mesh, thermo, U, phi));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fluidPhaseChangeModel::~fluidPhaseChangeModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


Foam::tmp<Foam::volScalarField>
Foam::fluidPhaseChangeModel::speciesSource(const word& specieName) const
{
    const label specieI = thermo_.composition().species().find(specieName);
    if (specieI != -1)
    {
        return speciesSource(specieI);
    }
    else
    {
        return zeroSpeciesSource_;
    }
}


bool Foam::fluidPhaseChangeModel::checkMassConservation
(
    const scalar relTol,
    const scalar absTol
) const
{
    if (!active_)
    {
        return true;
    }

    const speciesTable& species = thermo_.composition().species();
    volScalarField sumSpecies
    (
        IOobject
        (
            "sumSpeciesSource",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimMass/dimVolume/dimTime, Zero)
    );

    forAll(species, i)
    {
        sumSpecies += speciesSource(i);
    }

    tmp<volScalarField> tNetMass = massSource();
    const volScalarField& netMass = tNetMass();

    volScalarField diff = mag(sumSpecies - netMass);
    const scalar maxDiff = max(diff).value();

    const scalar maxSpeciesRate = max(mag(sumSpecies)).value();
    const scalar maxNetMassRate = max(mag(netMass)).value();
    const scalar maxRate = max(maxSpeciesRate, maxNetMassRate);

    // Relative discrepancy against peak mass transfer rate in domain
    const scalar relDiff = (maxRate > absTol) ? (maxDiff / maxRate) : maxDiff;

    if (relDiff > relTol && maxDiff > absTol)
    {
        #ifdef FULLDEBUG
        FatalErrorInFunction
            << "Mass conservation discrepancy in phaseChange model on region " << mesh_.name()
            << ": max|sum(speciesSource) - massSource| = " << maxDiff
            << " kg/(m^3 s) (relative error: " << relDiff
            << " > relTol " << relTol << ")"
            << exit(FatalError);
        #else
        WarningInFunction
            << "Mass conservation discrepancy in phaseChange model on region " << mesh_.name()
            << ": max|sum(speciesSource) - massSource| = " << maxDiff
            << " kg/(m^3 s) (relative error: " << relDiff
            << " > relTol " << relTol << ")" << endl;
        #endif
        return false;
    }

    return true;
}


// ************************************************************************* //
