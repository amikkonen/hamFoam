// Inner cells
forAll(mesh.cellZones(), cellZone) {
    // Set values based on temperature
    const labelList& selectedCells(mesh.cellZones()[cellZone]);

    

    // Internal cells   
	forAll(selectedCells, loopIndex) {
        const label& cellIndex = selectedCells[loopIndex];
        rho[cellIndex]      = rho_functions[cellZone].value(0);
        cp[cellIndex]       = cp_functions[cellZone].value(0);
    }
}


// Boundaries
volScalarField::Boundary& rhoBf     = rho.boundaryFieldRef();
volScalarField::Boundary& cpBf      = cp.boundaryFieldRef();

forAll(boundaryZoneMapping, patchIndex) 
{
    const polyPatch &ppatch = mesh.boundaryMesh()[patchIndex];
    // Only update if something to update.
    // Remember that createFields.H has a similar line.
    if (!isA<emptyPolyPatch>(ppatch) && !isA<wedgePolyPatch>(ppatch))
    {
        forAll(boundaryZoneMapping[patchIndex], patchFaceIndex)
        {
            const label& cellZoneIndex = boundaryZoneMapping[patchIndex][patchFaceIndex];
            rhoBf[patchIndex][patchFaceIndex] = rho_functions[cellZoneIndex].value(0);
            cpBf[patchIndex][patchFaceIndex]  = cp_functions[cellZoneIndex].value(0);
        }
    }
}

// Processor boundaries.
// TODO! Check if other boundaries are affected.
rho.correctBoundaryConditions();
cp.correctBoundaryConditions();
