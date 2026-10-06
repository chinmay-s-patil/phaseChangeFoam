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
    active_(dict_.getOrDefault<bool>("active", false))
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
            if (subDict.getOrDefault<bool>("active", false))
            {
                subDict.readIfPresent("phaseChangeMode", modelType);
                if (modelType == "none" || modelType == "off")
                {
                    modelType = "none";
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


// ************************************************************************* //
