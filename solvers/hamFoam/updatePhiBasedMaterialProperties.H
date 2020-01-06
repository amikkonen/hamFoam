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
        k[cellIndex]        = max(k_functions[cellZone].value(phi[cellIndex]),VSMALL);
        w[cellIndex]        = max(w_functions[cellZone].value(phi[cellIndex]),VSMALL);
        delta_p[cellIndex]  = max(delta_p_functions[cellZone].value(phi[cellIndex]),VSMALL);
        Dw[cellIndex]       = max(Dw_functions[cellZone].value(phi[cellIndex]),VSMALL);
        xi[cellIndex]       = max(xi_functions[cellZone].value(phi[cellIndex]),VSMALL); 
    }
}


// Boundaries
volScalarField::Boundary& phiBf     = phi.boundaryFieldRef();
volScalarField::Boundary& kBf       = k.boundaryFieldRef();
volScalarField::Boundary& wBf       = w.boundaryFieldRef();
volScalarField::Boundary& DwBf      = Dw.boundaryFieldRef();
volScalarField::Boundary& delta_pBf = delta_p.boundaryFieldRef();
volScalarField::Boundary& xiBf      = xi.boundaryFieldRef();

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
            kBf[patchIndex][patchFaceIndex] 
                = max(k_functions[cellZoneIndex].value(phiBf[patchIndex][patchFaceIndex]),VSMALL);
            wBf[patchIndex][patchFaceIndex] 
                = max(w_functions[cellZoneIndex].value(phiBf[patchIndex][patchFaceIndex]),VSMALL);
            delta_pBf[patchIndex][patchFaceIndex] 
                = max(delta_p_functions[cellZoneIndex].value(phiBf[patchIndex][patchFaceIndex]),VSMALL);
            DwBf[patchIndex][patchFaceIndex] 
                = max(Dw_functions[cellZoneIndex].value(phiBf[patchIndex][patchFaceIndex]),VSMALL); 
            xiBf[patchIndex][patchFaceIndex] 
                = max(xi_functions[cellZoneIndex].value(phiBf[patchIndex][patchFaceIndex]),VSMALL); 
        }
    }
}

// Processor boundaries.
// TODO! Check if other boundaries are affected.
k.correctBoundaryConditions();
w.correctBoundaryConditions();
delta_p.correctBoundaryConditions();
Dw.correctBoundaryConditions();
xi.correctBoundaryConditions();

