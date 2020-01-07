/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    hamFoam

Description
    Solves transient heat and moisture transport in porous media. 
    Diffusion and capillary transport. No air movement.
    
Author
    Antti Mikkonen, a.mikkonen@iki.fi, 2019


\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvOptions.H"
#include "bound.H"
#include "pimpleControl.H"
#include "wedgePolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"

    #include "createTime.H"
    #include "createMesh.H"

    pimpleControl pimple(mesh);

    #include "createFields.H"
    #include "setConstantMaterialProperties.H"
    #include "updateTBasedMaterialProperties.H"
    
//    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting calculation\n" << endl;

    while (runTime.run())
    {
        // TODO! Adjust timestep etc here...
        
        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;


        while (pimple.loop())
        {

            #include "updatePhiBasedMaterialProperties.H"

            fvScalarMatrix TEqn
            (
                fvm::ddt(rho*cp+cw*w,T) 
                - fvm::laplacian(k + hw*delta_p*phi*dpsDT, T)
                // separated
//                - fvm::laplacian(k, T)
//                - fvm::laplacian(hw*delta_p*phi*dpsDT, T)
                - fvc::laplacian(hw*delta_p*ps, phi)
             ==
                rho*cp*fvOptions(T)
            );
            TEqn.relax();
            fvOptions.constrain(TEqn);
            TEqn.solve();
            fvOptions.correct(T);
//            bound(T, TMin); //VSMALL, upped bound
            #include "updateTBasedMaterialProperties.H"

            fvScalarMatrix phiEqn
            (
                fvm::ddt(xi,phi) 
                - fvc::laplacian(delta_p*phi*dpsDT, T)
                - fvm::laplacian(Dw*xi + delta_p*ps, phi)

                 // Separated. Unstable?
//                - fvm::laplacian(delta_p*ps, phi)
//                - fvm::laplacian(Dw*xi, phi)
//             ==
//                xi*fvOptions(phi) // add suitable options
            );

            phiEqn.relax();
            fvOptions.constrain(phiEqn);
            phiEqn.solve();
            fvOptions.correct(phi);


            // Limit
//            phi = min(phi, phiMax);
//            phi = max(phi, phiMin);

        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;

}


// ************************************************************************* //














//            bound(phi, phiMin);//VSMALL, upped bound 1

            // Limits
//        	forAll(phi, index) {
//                T[index]   = min(T[index], TMax);
//                T[index]   = max(T[index], TMin);            
//                phi[index] = min(phi[index], phiMax);
//                phi[index] = max(phi[index], phiMin);
//            }
//            T   = min(T.primitiveField(), TMax.value());
////            T   = max(T, TMin);            
//            phi = min(phi, phiMax);
//            phi = max(phi, phiMin);

////            T.correctBoundaryConditions();
//            phi.correctBoundaryConditions();
//            phi.boundaryFieldRef() = min(phi.boundaryField(), phiMax);
//            phi.boundaryFieldRef() = max(phi.boundaryField(), phiMin);

