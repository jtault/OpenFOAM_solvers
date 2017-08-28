// Particle 1
volVectorField U_New = U + gammap1/Csalt*fvc::grad(Csalt);

surfaceScalarField phiNew = linearInterpolate(U_New) & mesh.Sf();

tmp<fvScalarMatrix> Cpart1Eqn
(
 fvm::ddt(Cpart1)
 +fvm::div(phiNew,Cpart1)
 -fvm::laplacian(Dpart,Cpart1)
 );

solve(Cpart1Eqn());

// Particle 2
volVectorField U_New = U + gammap2/Csalt*fvc::grad(Csalt);

surfaceScalarField phiNew = linearInterpolate(U_New) & mesh.Sf();

tmp<fvScalarMatrix> Cpart2Eqn
(
 fvm::ddt(Cpart2)
 +fvm::div(phiNew,Cpart2)
 -fvm::laplacian(Dpart,Cpart2)
 );

solve(Cpart2Eqn());

// Particle 3
volVectorField U_New = U + gammap3/Csalt*fvc::grad(Csalt);

surfaceScalarField phiNew = linearInterpolate(U_New) & mesh.Sf();

tmp<fvScalarMatrix> Cpart3Eqn
(
 fvm::ddt(Cpart3)
 +fvm::div(phiNew,Cpart3)
 -fvm::laplacian(Dpart,Cpart3)
 );

solve(Cpart3Eqn());

// Particle 4
volVectorField U_New = U + gammap4/Csalt*fvc::grad(Csalt);

surfaceScalarField phiNew = linearInterpolate(U_New) & mesh.Sf();

tmp<fvScalarMatrix> Cpart4Eqn
(
 fvm::ddt(Cpart4)
 +fvm::div(phiNew,Cpart4)
 -fvm::laplacian(Dpart,Cpart4)
 );

solve(Cpart4Eqn());

// Particle 5
volVectorField U_New = U + gammap5/Csalt*fvc::grad(Csalt);

surfaceScalarField phiNew = linearInterpolate(U_New) & mesh.Sf();

tmp<fvScalarMatrix> Cpart5Eqn
(
 fvm::ddt(Cpart5)
 +fvm::div(phiNew,Cpart5)
 -fvm::laplacian(Dpart,Cpart5)
 );

solve(Cpart5Eqn());
