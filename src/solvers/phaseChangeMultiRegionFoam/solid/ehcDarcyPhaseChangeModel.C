/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Website:  www.openfoam.com
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

#include "ehcDarcyPhaseChangeModel.H"
#include "calculatedFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(ehcDarcyPhaseChangeModel, 0);
    addToRunTimeSelectionTable(phaseChangeModel, ehcDarcyPhaseChangeModel, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ehcDarcyPhaseChangeModel::ehcDarcyPhaseChangeModel
(
    const fvMesh& mesh,
    const basicThermo& thermo
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


void Foam::ehcDarcyPhaseChangeModel::readPorosityDict()
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
                for (const word& k : porDict.toc())
                {
                    if (k != "Cu" && k != "q" && k != "A_cu" && k != "eps")
                    {
                        FatalIOErrorInFunction(porDict)
                            << "Unknown key '" << k << "' in porosity dict."
                            << exit(FatalIOError);
                    }
                }
                Cu_ = porDict.getOrDefault<scalar>("Cu", porDict.lookupOrDefault<scalar>("A_cu", 1.0e5));
                q_  = porDict.getOrDefault<scalar>("q", porDict.lookupOrDefault<scalar>("eps", 0.001));
            }

            if (pcDict.found("buoyancy"))
            {
                const dictionary& buoyDict = pcDict.subDict("buoyancy");
                for (const word& k : buoyDict.toc())
                {
                    if (k != "beta")
                    {
                        FatalIOErrorInFunction(buoyDict)
                            << "Unknown key '" << k << "' in buoyancy dict."
                            << exit(FatalIOError);
                    }
                }
                beta_ = buoyDict.lookupOrDefault<scalar>("beta", 0.0);
            }
        }
    }
}


Foam::tmp<Foam::volScalarField> Foam::ehcDarcyPhaseChangeModel::momentumSp() const
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
    const scalarField& alphaCells = phaseFraction_.primitiveField();
    scalarField& coeffCells = coeff.primitiveFieldRef();

    forAll(coeffCells, celli)
    {
        scalar a = alphaCells[celli];
        // Darcy mushy zone drag coefficient: Cu * (1 - alpha)^2 / (alpha^3 + q)
        coeffCells[celli] = Cu_ * sqr(1.0 - a) / (pow3(a) + q_);
    }

    return tCoeff;
}

// ************************************************************************* //
