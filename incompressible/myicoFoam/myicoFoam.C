/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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
    icoFoam

Description
    Transient solver for incompressible, laminar flow of Newtonian fluids.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "pisoControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline scalar ip(const volVectorField& u1, const volVectorField& u2){
  return (u1 & u2)().weightedAverage(u1.mesh().V()).value();
}

inline scalar norm(const volVectorField& u) {
  return ::sqrt(ip(u,u));
}

inline scalar scalar_norm(const volScalarField& u) {
  return ::sqrt((sqr(u))().weightedAverage(u.mesh().V()).value());
}

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    pisoControl piso(mesh);

    #include "createFields.H"
    #include "initContinuityErrs.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "CourantNo.H"

        // Momentum predictor

        fvVectorMatrix UEqn
        (
            fvm::ddt(U)
          + fvm::div(phi, U)
          - fvm::laplacian(nu, U)
        );

	volVectorField Old_U = U;

        if (piso.momentumPredictor())
        {
            solve(UEqn == -fvc::grad(p));
        }

        // --- PISO loop
        while (piso.correct())
        {
            volScalarField rAU(1.0/UEqn.A());

            volVectorField HbyA("HbyA", U);
            HbyA = rAU*UEqn.H();
            surfaceScalarField phiHbyA
            (
                "phiHbyA",
                (fvc::interpolate(HbyA) & mesh.Sf())
              + fvc::interpolate(rAU)*fvc::ddtCorr(U, phi)
            );

            adjustPhi(phiHbyA, U, p);

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

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
	double dU_dt = norm(U - Old_U) / runTime.deltaT().value();
	Info<< "norm(dU_dt) = " << dU_dt << endl;

	volVectorField new_U = fvc::reconstruct(phi);
	new_U.correctBoundaryConditions();

	Info << "norm(new_U) = " << norm(new_U) << endl;
	Info << "norm(U) = " << norm(U) << endl;
	Info << "norm(U - new_U) = " << norm(U - new_U) << endl;
	//Info << "Interpolated velocity error = "
	//     << (sqrt(sum(sqr(fvc::flux(new_U) - phi)))/sum(mesh.magSf())).value()
	//     << endl;  

	Info << "norm(div(U)) = " << scalar_norm(fvc::div(U)) << endl;
	Info << "norm(div(phi) = " << scalar_norm(fvc::div(phi)) << endl;
        Info << "max(div(U)) = " << max(fvc::div(U)) << endl;
	Info << "max(div(phi)) = " << max(fvc::div(phi)) << endl;
	//Info << "max(phi) = " << max(phi) << endl;
	//Info << "max(phi_from_U) = " << max(phi_from_U) << endl;
	//Info << "max(phi - phi_from_U) = " << max(phi - phi_from_U) << endl;
	//Info << "max(fvc::div(fvc::interpolate(U))) = " << max(fvc::div(fvc::interpolate(U))) << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
