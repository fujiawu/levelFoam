    
    Info<< "Reading constants\n" << endl;
    IOdictionary globalConstants
    (
        IOobject
        (
            "globalConstants",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    Info<< "Reading kinematic viscosity nu\n" << endl;
    dimensionedScalar nu
    (
        globalConstants.lookup("nu")
    );

    Info<< "Reading Schmidt number sc\n" << endl;
    dimensionedScalar sc
    (

        globalConstants.lookup("sc")
    );

    Info<< "Reading Lewis number le\n" << endl;
    dimensionedScalar le
    (
        globalConstants.lookup("le")
    );
     
    Info<< "Reading thermal expansion ratio theta\n" << endl;
    dimensionedScalar theta
    (
        globalConstants.lookup("theta")
    );

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

    Info<< "Reading field c\n" << endl;
    volScalarField c
    (
        IOobject
        (
            "c",
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

   Info<< "Creating field psi\n" << endl;
    volScalarField psi
    (
       IOobject
       (
           "psi",
           runTime.timeName(),
           mesh,
           IOobject::NO_READ,
           IOobject::NO_WRITE
       ),
       Psi(c,theta)
    );

    Info<< "Creating field rho\n" << endl;
    volScalarField rho
    (
       IOobject
       (
           "rho",
           runTime.timeName(),
           mesh,
           IOobject::NO_READ,
           IOobject::AUTO_WRITE
       ),
       Rho(p,psi)
    );

 
