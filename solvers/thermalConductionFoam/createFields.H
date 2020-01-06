Info<< "Reading field T\n" << endl;

volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading transportProperties\n" << endl;

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);


////////////////////////////////////////////////////////////////////////
// Read cellZones and check number of cells
////////////////////////////////////////////////////////////////////////

scalar totalCellsInZones = 0;

forAll(mesh.cellZones(), cellZone) {
    word zoneName = mesh.cellZones()[cellZone].name();
    scalar ncells = returnReduce(mesh.cellZones()[cellZone].size(), sumOp<label>());
    totalCellsInZones += ncells;
    Info << "Found zone " << zoneName << " with " << ncells << " cells." << endl;
}
// Not needed because ncells reduced? 
//reduce(totalCellsInZones, sumOp<scalar>());

// Check that cells in zones match total cell number
// NOTE! Assumes zones are not overlapping
if (totalCellsInZones == mesh.globalData().nTotalCells()) { 
    Info << "Total cells in zones matches total cells in mesh: " 
         << totalCellsInZones << endl;
} else {
    FatalErrorInFunction
                << "Total cells in zones: " << totalCellsInZones 
                << " does not match mesh size "
                << mesh.globalData().nTotalCells()
                << exit(FatalError);
}

////////////////////////////////////////////////////////////////////////
// Create material properties, values in updateMaterialProperties.H
////////////////////////////////////////////////////////////////////////

volScalarField rho
(
     IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimDensity, 0)
);

volScalarField cp
(
     IOobject
    (
        "cp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimEnergy/dimMass/dimTemperature, 0)
);

volScalarField k
(
     IOobject
    (
        "k",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimPower/dimTemperature/dimLength, 0)
);

////////////////////////////////////////////////////////////////////////
// Init material functions
////////////////////////////////////////////////////////////////////////

PtrList<Function1<scalar>> rho_functions;
PtrList<Function1<scalar>> cp_functions;
PtrList<Function1<scalar>> k_functions;

forAll(mesh.cellZones(), cellZone) {
    const word& materialName = mesh.cellZones()[cellZone].name();

    const dictionary& materialDict = transportProperties.subDict(materialName);

    // density (kg/m^3)
    autoPtr<Function1<scalar>> rho_function1
    (
        Function1<scalar>::New
        (
            "rho",
            materialDict
        )
    );
    rho_functions.append(rho_function1);

    // specific heat capacity (J/(kgK)) 
    autoPtr<Function1<scalar>> cp_function1
    (
        Function1<scalar>::New
        (
            "cp",
            materialDict
        )
    );
    cp_functions.append(cp_function1);

    // thermal conductivity (W/(mK))
    autoPtr<Function1<scalar>> k_function1
    (
        Function1<scalar>::New
        (
            "k",
            materialDict
        )
    );
    k_functions.append(k_function1);

    // Info
    Info << materialName << endl
         << "\trho type " << rho_functions[cellZone].type() << endl
         << "\tcp  type " << cp_functions[cellZone].type() << endl
         << "\tk   type " << k_functions[cellZone].type() << endl;
}

////////////////////////////////////////////////////////////////////////
// Create mapping for boundary materials. 
// This is a ugly way of doing this...
// Hopefulle there is something better available.
////////////////////////////////////////////////////////////////////////

// Contains the cellZone index of the faceCells of patches. 
// Example usage:
// const label& cellZoneIndex = boundaryZoneMapping[patchIndex][patchFaceIndex];
DynamicList<labelList> boundaryZoneMapping(0);

forAll(mesh.boundaryMesh(), patchi)
{
    labelList internalFieldNextToPatch = mesh.boundaryMesh()[patchi].faceCells(); 
    labelList patchMapping(internalFieldNextToPatch.size());

    const polyPatch &ppatch = mesh.boundaryMesh()[patchi];
    // Only add to mapping if something to map.
    // Remember that updateMaterialProperties.H has a similar line.
    if (!isA<emptyPolyPatch>(ppatch) && !isA<wedgePolyPatch>(ppatch))
    {
        forAll(internalFieldNextToPatch, internalFieldNextToPatchk)
        {
            const label& internalCellIndex = internalFieldNextToPatch[internalFieldNextToPatchk];
            // Check if index in zone.
            bool breakFlag = false;
            forAll(mesh.cellZones(), cellZone) {
                const labelList& selectedCells(mesh.cellZones()[cellZone]);
                // Internal cells   
	            forAll(selectedCells, loopIndex) {
                    const label& cellIndexInZone = selectedCells[loopIndex];
                    // Check if match.
                    // If match, store index and break out of two loops.
                    if (internalCellIndex == cellIndexInZone) 
                    {
                        patchMapping[internalFieldNextToPatchk] = cellZone;
                        breakFlag = true;
                        break;
                    }
                }
                if (breakFlag) {break;}
            }
        }
    }

    // One patch mapping ready. Store.
    boundaryZoneMapping.append(patchMapping);    
}

//////////////////////////////////////
// Tesing of boundaryZoneMapping
//////////////////////////////////////
//forAll(boundaryZoneMapping, patchIndex) 
//{
//    //Info << boundaryZoneMapping[patchIndex] << endl;
//    Info << patchIndex << endl;  
//    forAll(boundaryZoneMapping[patchIndex], patchFaceIndex)
//    {
//        const label& cellZoneIndex = boundaryZoneMapping[patchIndex][patchFaceIndex];
//        Info << "    " << cellZoneIndex << endl;    
//    }
//}


#include "createFvOptions.H"

