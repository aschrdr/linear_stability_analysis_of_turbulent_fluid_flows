/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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
    icoFoam

Group
    grpIncompressibleSolvers

Description
    Transient solver for incompressible, laminar flow of Newtonian fluids.

    \heading Solver details
    The solver uses the PISO algorithm to solve the continuity equation:

        \f[
            \div \vec{U} = 0
        \f]

    and momentum equation:

        \f[
            \ddt{\vec{U}}
          + \div \left( \vec{U} \vec{U} \right)
          - \div \left(\nu \grad \vec{U} \right)
          = - \grad p
        \f]

    Where:
    \vartable
        \vec{U} | Velocity
        p       | Pressure
    \endvartable

    \heading Required fields
    \plaintable
        U       | Velocity [m/s]
        p       | Kinematic pressure, p/rho [m2/s2]
    \endplaintable

\*---------------------------------------------------------------------------*/
#include <iostream>
#include <vector>

#include "fvCFD.H"
#include "pisoControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

std::vector<double> linspace(double start, double end, int num) {
    std::vector<double> result;
    double step = (end - start) / (num - 1);

    for (int i = 0; i < num; ++i) {
        result.push_back(start + i * step);
    }

    return result;
} // Sourced from the internet: how to get numpy linspace equivalent in c++

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Transient solver for incompressible, laminar flow"
        " of Newtonian fluids."
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

    pisoControl piso(mesh);

    #include "createFields.H"
    #include "initContinuityErrs.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    volVectorField v
    (
     IOobject
     (
      "U",
      runTime.timeName(),
      mesh,
      IOobject::MUST_READ,
      IOobject::AUTO_WRITE
     ),
     mesh
    );

    int n = 10;
    int j = 1;
    std::vector<double> timeSteps = linspace(0,0.001,n+1);
    SquareMatrix<scalar> Hessenberg(n+1,0);

    double vmag;
    vmag = Foam::sqrt(sum(magSqr(v)).value());
    v = v/vmag;

    std::cout << Foam::sqrt(sum(magSqr(v)).value()) << std::endl;

    runTime.write();

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "CourantNo.H"

        // Momentum predictor

        fvVectorMatrix UEqn
        (
            fvm::ddt(U)
          + fvm::div(phiB, U)
	  + fvc::div(phi, UB)
          - fvm::laplacian(nu, U)
        );

        if (piso.momentumPredictor())
        {
            solve(UEqn == -fvc::grad(p));
        }

        // --- PISO loop
        while (piso.correct())
        {
            volScalarField rAU(1.0/UEqn.A());
            volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
            surfaceScalarField phiHbyA
            (
                "phiHbyA",
                fvc::flux(HbyA)
              + fvc::interpolate(rAU)*fvc::ddtCorr(U, phi)
            );

            adjustPhi(phiHbyA, U, p);

            // Update the pressure BCs to ensure flux consistency
            constrainPressure(p, U, phiHbyA, rAU);

            // Non-orthogonal pressure corrector loop
            while (piso.correctNonOrthogonal())
            {
                // Pressure corrector

                fvScalarMatrix pEqn
                (
                    fvm::laplacian(rAU, p) == fvc::div(phiHbyA)
                );

                pEqn.setReference(pRefCell, pRefValue);

                pEqn.solve(mesh.solver(p.select(piso.finalInnerIter())));

                if (piso.finalNonOrthogonalIter())
                {
                    phi = phiHbyA - pEqn.flux();
                }
            }

            #include "continuityErrs.H"

            U = HbyA - rAU*fvc::grad(p);
            U.correctBoundaryConditions();
        }

	// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	
	double dt = static_cast<double>(runTime.value());

	if (std::abs(timeSteps[j] - dt) < 0.0000001) {
		for (int i=0; i < j; i++){
			volVectorField vi
			(
			 IOobject
			 (
			  "U",
			  runTime.timeName(timeSteps[i]),
			  mesh,
			  IOobject::MUST_READ,
			  IOobject::NO_WRITE
			 ),
			 mesh
			);

			// std::cout << Foam::sqrt(sum(magSqr(vi)).value()) << std::endl;

			volVectorField prodField = cmptMultiply(U, vi);
			Hessenberg[i][j-1] = sum(prodField.component(0) + prodField.component(1) + prodField.component(2)).value();
		}

		for (int i=0; i < j; i++){
			volVectorField vi
			(
			 IOobject
			 (
			  "U",
			  runTime.timeName(timeSteps[i]),
			  mesh,
			  IOobject::MUST_READ,
			  IOobject::NO_WRITE
			 ),
			 mesh
			);

			U = U - Hessenberg[i][j-1]*vi;
		}

		Hessenberg[j][j-1] = Foam::sqrt(sum(magSqr(U))).value();
		U = U/Hessenberg[j][j-1];

		// std::cout << Foam::sqrt(sum(magSqr(U)).value()) << std::endl;

		j += 1;
	}

	// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    // Write the Hessenberg
    
    #include "OFstream.H"

    fileName name = ("Hessenberg.csv");
    OFstream OS(name);

    for(int i=0; i<n+1; i++){
	    for(int j=0; j<n+1; j++){
		    if(OS.opened() && Pstream::master()){
			    OS << Hessenberg[i][j];
			    if(j != n){
				    OS << ",";
			    }
		    }
	    }
	    OS << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
