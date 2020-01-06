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
    steadyHamFoam

Description
    Solves steady state heat and moisture transport in porous media. 
    Diffusion and capillary transport. No air movement.
    
Author
    Antti Mikkonen, a.mikkonen@iki.fi, 2019


\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvOptions.H"
#include "bound.H"
#include "simpleControl.H"

#include "wedgePolyPatch.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"

    #include "createTime.H"
    #include "createMesh.H"

    simpleControl simple(mesh);

    #include "createFields.H"
    #include "setConstantMaterialProperties.H"
    #include "updateTBasedMaterialProperties.H"
    
//    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nnStarting calculation\n" << endl;


    while (simple.loop(runTime))
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        while (simple.correctNonOrthogonal())
        {

            #include "updatePhiBasedMaterialProperties.H"

            fvScalarMatrix TEqn
            (
                - fvm::laplacian(k, T)
                - fvm::laplacian(hw*delta_p*phi*dpsDT, T)
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
                - fvc::laplacian(delta_p*phi*dpsDT, T)
                - fvm::laplacian(delta_p*ps, phi)
                - fvm::laplacian(Dw*xi, phi) // combine with abowe?
//             ==
//                xi*fvOptions(phi) // add suitable options
            );

            phiEqn.relax();
            fvOptions.constrain(phiEqn);
            phiEqn.solve();
            fvOptions.correct(phi);


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


