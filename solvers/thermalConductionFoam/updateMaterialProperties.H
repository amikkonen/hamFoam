// Inner cells
forAll(mesh.cellZones(), cellZone) {
    // Set values based on temperature
    const labelList& selectedCells(mesh.cellZones()[cellZone]);

    // TODO! 
    // MOD harmonic scheme for zeros.
    // NOTE! SMALL is too large

    // Internal cells   
	forAll(selectedCells, loopIndex) {
        const label& cellIndex = selectedCells[loopIndex];
        rho[cellIndex] = max(rho_functions[cellZone].value(T[cellIndex]),VSMALL);
        cp[cellIndex]  = max(cp_functions[cellZone].value(T[cellIndex]),VSMALL);        
        k[cellIndex]   = max(k_functions[cellZone].value(T[cellIndex]),VSMALL);
    }
}

// Boundaries
volScalarField::Boundary& TBf   = T.boundaryFieldRef();
volScalarField::Boundary& rhoBf = rho.boundaryFieldRef();
volScalarField::Boundary& cpBf  = cp.boundaryFieldRef();
volScalarField::Boundary& kBf   = k.boundaryFieldRef();

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
            rhoBf[patchIndex][patchFaceIndex] 
                = max(rho_functions[cellZoneIndex].value(TBf[patchIndex][patchFaceIndex]),VSMALL);
            cpBf[patchIndex][patchFaceIndex] 
                = max(cp_functions[cellZoneIndex].value(TBf[patchIndex][patchFaceIndex]),VSMALL);
            kBf[patchIndex][patchFaceIndex] 
                = max(k_functions[cellZoneIndex].value(TBf[patchIndex][patchFaceIndex]),VSMALL);
        }
    }
}

// Processor boundaries.
rho.correctBoundaryConditions();
cp.correctBoundaryConditions();
k.correctBoundaryConditions();


