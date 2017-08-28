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
    sangwooSimpleFoam

Description
    Steady-state solver for incompressible, turbulent flow

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "simpleControl.H"
#include "fvOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline scalar ip(const volVectorField& u1, const volVectorField& u2) {
  return (u1 & u2)().weightedAverage(u1.mesh().V()).value();
}

inline scalar norm(const volVectorField& u) {
  return ::sqrt(ip(u,u));
}

int main(int argc, char *argv[])
{
    #include "postProcess.H"
  
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "createFvOptions.H"
    #include "initContinuityErrs.H"

    turbulence->validate();
  
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    int myCounter = 1;

    while (runTime.run())
      {
        // Update the concentration numSteps times
        for (int n = 0; n < numSteps.value(); ++n) {
          runTime++;
          Info<< "Time = " << runTime.timeName() << nl << endl;

          #include "CsaltEqn.H"

          runTime.write();

          Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
              << "  ClockTime = " << runTime.elapsedClockTime() << " s"
              << nl << endl;
        }

        // SOLVE THE FLOW FIELD
        while (true)
          {
            volVectorField Old_U = U;
            // --- Pressure-velocity SIMPLE corrector
            #include "UEqn.H"
            #include "pEqn.H"
	    laminarTransport.correct();
            turbulence->correct();

            double dt=runTime.deltaT().value();
            double dU_dt = norm(U - Old_U) / dt;
            Info<< "norm(dU_dt) = " << dU_dt << endl;
            Info<< "counter = " << myCounter << endl;
            if ((dU_dt<myTol.value() && myCounter>minIters.value()) || myCounter>numIters.value())
              {
                break;
              }
            ++myCounter;
          }

        myCounter = 1;
      }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
