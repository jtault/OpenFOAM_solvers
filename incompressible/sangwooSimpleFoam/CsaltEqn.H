while (simple.correctNonOrthogonal())
{
    fvScalarMatrix CsaltEqn
    (
       fvm::ddt(Csalt)
       +fvm::div(phi,Csalt)
       -fvm::laplacian(Dsalt,Csalt)
    );

    CsaltEqn.relax();
    fvOptions.constrain(CsaltEqn);
    CsaltEqn.solve();
    fvOptions.correct(Csalt);
}
