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
    const solidThermo& thermo
)
:
    phaseChangeModel(mesh, thermo, false),
    beta_(0.0),
    Cu_(1.0e5),
    q_(0.001)
{

    readPorosityDict();

    // Validate patch types
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
    checkCalculatedPatches(Cp_);
    checkCalculatedPatches(rho_);
    checkCalculatedPatches(k_);
}


void Foam::enthalpyPorosityPhaseChangeModel::readPorosityDict()
{
    IOobject dictIO
    (
        "phaseChangeDict",
        mesh_.time().constant(),
        mesh_,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    );

    if (dictIO.typeHeaderOk<dictionary>(true))
    {
        IOdictionary phaseChangeDict(dictIO);
        if (phaseChangeDict.found("phaseChange"))
        {
            const dictionary& pcDict = phaseChangeDict.subDict("phaseChange");

            // Allow optional convection dict override
            if (pcDict.found("convection"))
            {
                suppressConvection_ =
                    pcDict.subDict("convection").getOrDefault<bool>("suppress", false);
            }

            // Read momentum damping parameters
            if (pcDict.found("porosity"))
            {
                const dictionary& porDict = pcDict.subDict("porosity");
                Cu_ = porDict.lookupOrDefault<scalar>("Cu", 1.0e5);
                q_  = porDict.lookupOrDefault<scalar>("q", 0.001);
            }

            if (pcDict.found("buoyancy"))
            {
                const dictionary& buoyDict = pcDict.subDict("buoyancy");
                beta_ = buoyDict.lookupOrDefault<scalar>("beta", 0.0);
            }
        }
    }
}


Foam::tmp<Foam::volVectorField> Foam::enthalpyPorosityPhaseChangeModel::momentumSource() const
{
    auto tSource = tmp<volVectorField>::New
    (
        IOobject
        (
            "momentumSource",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector("zero", dimForce/dimVolume, Zero)
    );

    if (!active_ || Cu_ <= 0.0)
    {
        return tSource;
    }

    // Lookup velocity field U if present in region
    if (mesh_.foundObject<volVectorField>("U"))
    {
        const volVectorField& U = mesh_.lookupObject<volVectorField>("U");
        volVectorField& source = tSource.ref();

        const scalarField& alphaCells = phaseFraction_.primitiveField();
        const vectorField& Ucells = U.primitiveField();
        vectorField& sourceCells = source.primitiveFieldRef();

        forAll(sourceCells, celli)
        {
            scalar a = alphaCells[celli];
            // Darcy mushy zone drag: - Cu * (1 - alpha)^2 / (alpha^3 + q) * U
            scalar A_mush = Cu_ * sqr(1.0 - a) / (pow3(a) + q_);
            sourceCells[celli] = - A_mush * Ucells[celli];
        }
    }

    return tSource;
}

// ************************************************************************* //
