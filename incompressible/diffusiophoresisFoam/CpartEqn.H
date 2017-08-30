while (simple.correctNonOrthogonal())
{

    // Particle 1
    volVectorField U_New1 = U + gammap1/Csalt*fvc::grad(Csalt);
    surfaceScalarField phiNew1 = linearInterpolate(U_New1) & mesh.Sf();

    fvScalarMatrix Cpart1Eqn
    (
       fvm::ddt(Cpart1)
       +fvm::div(phiNew1,Cpart1)
       -fvm::laplacian(Dpart1,Cpart1)
    );

    Cpart1Eqn.relax();
    fvOptions.constrain(Cpart1Eqn);
    Cpart1Eqn.solve();
    fvOptions.correct(Cpart1);

    // Particle 2
    volVectorField U_New2 = U + gammap2/Csalt*fvc::grad(Csalt);
    surfaceScalarField phiNew2 = linearInterpolate(U_New2) & mesh.Sf();

    fvScalarMatrix Cpart2Eqn
    (
       fvm::ddt(Cpart2)
       +fvm::div(phiNew2,Cpart2)
       -fvm::laplacian(Dpart2,Cpart2)
    );

    Cpart2Eqn.relax();
    fvOptions.constrain(Cpart2Eqn);
    Cpart2Eqn.solve();
    fvOptions.correct(Cpart2);

    // Particle 3
    volVectorField U_New3 = U + gammap3/Csalt*fvc::grad(Csalt);
    surfaceScalarField phiNew3 = linearInterpolate(U_New3) & mesh.Sf();

    fvScalarMatrix Cpart3Eqn
    (
       fvm::ddt(Cpart3)
       +fvm::div(phiNew3,Cpart3)
       -fvm::laplacian(Dpart3,Cpart3)
    );

    Cpart3Eqn.relax();
    fvOptions.constrain(Cpart3Eqn);
    Cpart3Eqn.solve();
    fvOptions.correct(Cpart3);

    // Particle 4
    volVectorField U_New4 = U + gammap4/Csalt*fvc::grad(Csalt);
    surfaceScalarField phiNew4 = linearInterpolate(U_New4) & mesh.Sf();

    fvScalarMatrix Cpart4Eqn
    (
       fvm::ddt(Cpart4)
       +fvm::div(phiNew4,Cpart4)
       -fvm::laplacian(Dpart4,Cpart4)
    );

    Cpart4Eqn.relax();
    fvOptions.constrain(Cpart4Eqn);
    Cpart4Eqn.solve();
    fvOptions.correct(Cpart4);

    // Particle 5
    volVectorField U_New5 = U + gammap5/Csalt*fvc::grad(Csalt);
    surfaceScalarField phiNew5 = linearInterpolate(U_New5) & mesh.Sf();

    fvScalarMatrix Cpart5Eqn
    (
       fvm::ddt(Cpart5)
       +fvm::div(phiNew5,Cpart5)
       -fvm::laplacian(Dpart5,Cpart5)
    );

    Cpart5Eqn.relax();
    fvOptions.constrain(Cpart5Eqn);
    Cpart5Eqn.solve();
    fvOptions.correct(Cpart5);

}
