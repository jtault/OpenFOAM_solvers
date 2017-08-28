Info<< "Reading transportProperties\n" << endl;
IOdictionary transportProperties
(
 IOobject
 (
  "transportProperties",
  runTime.constant(),
  mesh,
  IOobject::MUST_READ,
  IOobject::NO_WRITE
  )
 );

dimensionedScalar nu(transportProperties.lookup("nu"));

Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
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

Info<< "Reading field Csalt\n" << endl;
volScalarField Csalt
(
    IOobject
    (
     "Csalt",
     runTime.timeName(),
     mesh,
     IOobject::MUST_READ,
     IOobject::AUTO_WRITE
     ),
    mesh
 );

Info<< "Reading field CsaltGrad\n" << endl;
volVectorField CsaltGrad
(
    IOobject
    (
     "CsaltGrad",
     runTime.timeName(),
     mesh,
     IOobject::MUST_READ,
     IOobject::AUTO_WRITE
     ),
    mesh
 );

Info<< "Reading field Cpart1\n" << endl;
volScalarField Cpart1
(
    IOobject
    (
     "Cpart1",
     runTime.timeName(),
     mesh,
     IOobject::MUST_READ,
     IOobject::AUTO_WRITE
     ),
    mesh
 );

Info<< "Reading field Cpart2\n" << endl;
volScalarField Cpart2
(
    IOobject
    (
     "Cpart2",
     runTime.timeName(),
     mesh,
     IOobject::MUST_READ,
     IOobject::AUTO_WRITE
     ),
    mesh
 );

Info<< "Reading field Cpart3\n" << endl;
volScalarField Cpart3
(
    IOobject
    (
     "Cpart3",
     runTime.timeName(),
     mesh,
     IOobject::MUST_READ,
     IOobject::AUTO_WRITE
     ),
    mesh
 );

Info<< "Reading field Cpart4\n" << endl;
volScalarField Cpart4
(
    IOobject
    (
     "Cpart4",
     runTime.timeName(),
     mesh,
     IOobject::MUST_READ,
     IOobject::AUTO_WRITE
     ),
    mesh
 );

Info<< "Reading field Cpart5\n" << endl;
volScalarField Cpart5
(
    IOobject
    (
     "Cpart5",
     runTime.timeName(),
     mesh,
     IOobject::MUST_READ,
     IOobject::AUTO_WRITE
     ),
    mesh
 );

dimensionedScalar Dsalt(transportProperties.lookup("Dsalt"));

dimensionedScalar Dpart1(transportProperties.lookup("Dpart1"));
dimensionedScalar gammap1(transportProperties.lookup("gammap1"));
dimensionedScalar Dpart2(transportProperties.lookup("Dpart2"));
dimensionedScalar gammap2(transportProperties.lookup("gammap2"));
dimensionedScalar Dpart3(transportProperties.lookup("Dpart3"));
dimensionedScalar gammap3(transportProperties.lookup("gammap3"));
dimensionedScalar Dpart4(transportProperties.lookup("Dpart4"));
dimensionedScalar gammap4(transportProperties.lookup("gammap4"));
dimensionedScalar Dpart5(transportProperties.lookup("Dpart5"));
dimensionedScalar gammap5(transportProperties.lookup("gammap5"));

dimensionedScalar timeFactor(transportProperties.lookup("timeFactor"));
dimensionedScalar myTol(transportProperties.lookup("myTol"));
dimensionedScalar numSteps(transportProperties.lookup("numSteps"));
dimensionedScalar numIters(transportProperties.lookup("numIters"));
dimensionedScalar minIters(transportProperties.lookup("minIters"));

#include "createPhi.H"

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, simple.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());

singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);

#include "createMRF.H"
