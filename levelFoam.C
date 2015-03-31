/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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
    levelFoam

Description

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
//#include "multivariateScheme.H"
//#include "ODESolvers.H"
//#include "fvIOoptionList.H"
//#include "meshTools.H"

// State Equation for Density
tmp<volScalarField> Rho(volScalarField& p, volScalarField& psi)
{
  return p * psi;

}

// State Equation for psi
tmp<volScalarField> Psi(volScalarField& c, dimensionedScalar& theta)
{

 // a constant to get the absolute density right
  dimensionedScalar alpha("alpha", 
      dimensionSet(0, -2, 2, 0, 0), 1);

 return alpha / ( 1 + c * (theta - 1));

}

// Reaction Rate Equation
tmp<volScalarField> Reaction(volScalarField& c)
{
  double A = 0.1;
  return A * c * (1-c);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:
int main(int argc, char *argv[])
{
  #include "setRootCase.H"
  #include "createTime.H"
  #include "createMesh.H"
  #include "createFields.H"
  #include "compressibleCreatePhi.H"
  #include "initContinuityErrs.H"
  #include "readTimeControls.H"
  #include "compressibleCourantNo.H"
  #include "setInitialDeltaT.H"
  //#include "createFvOptions.H"

  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
  while (runTime.loop())
  {

    Info<< "Time = " << runTime.timeName() << nl << endl;
 
    #include "readTimeControls.H"
    #include "readPISOControls.H"
    #include "compressibleCourantNo.H"
    #include "setDeltaT.H"
    //Info << "nOuterCorr:" << nOuterCorr << endl;    
    //Info << "nCorr:" << nCorr << endl;    
      
    // Reaction Eqn
 
    // Density Eqn
    fvScalarMatrix rhoEqn(fvm::ddt(rho) + fvc::div(phi));
    rhoEqn.solve();

    for (label ocorr=1; ocorr <= nOuterCorr; ocorr++)
    {
 
       // momentum Eqn
       fvVectorMatrix UEqn(fvm::ddt(rho, U) + fvm::div(phi, U)
        - fvm::laplacian(rho*nu, U));
       //UEqn.relax();
       solve(UEqn == - fvc::grad(p));

       // transport Eqn
       fvScalarMatrix cEqn(fvm::ddt(rho, c) + fvm::div(phi, c)
          - fvm::laplacian(rho*nu/sc, c));
       solve(cEqn);
   
       // Reaction Eqn, runTime.value();
       c = c + Reaction(c) * runTime.deltaT().value();
       psi = Psi(c, theta);


       volScalarField rAU(1.0/UEqn.A());
       volVectorField HbyA("HbyA", U);
       for (int corr=1; corr<=nCorr; corr++) 
       {
          // State Eqn
          rho = Rho(p, psi);
          HbyA = rAU*UEqn.H();
          surfaceScalarField phiHbyA("phiHbyA",
                  (fvc::interpolate(rho*HbyA) & mesh.Sf())
                  + fvc::interpolate(rho*rAU) * fvc::ddtCorr(rho, U, phi));

          adjustPhi(phiHbyA, U, p);    

          for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++) 
          {

             // Pressure Eqn
             fvScalarMatrix pEqn(fvm::laplacian(rho*rAU, p) 
                               == fvc::div(phiHbyA) + fvm::ddt(psi, p));
             pEqn.solve();

             if (nonOrth == nNonOrthCorr) phi = phiHbyA - pEqn.flux();
          }

          // Density Eqn
          fvScalarMatrix rhoEqn(fvm::ddt(rho) + fvc::div(phi));
          rhoEqn.solve();
          
          //#include "compressibleContinuityErrs.H"
       }

       U = HbyA - rAU*fvc::grad(p);
       U.correctBoundaryConditions();
       //DpDt = fvc::DDt(surfaceScalarField("phiU", phi/fvc::interpolate(rho)), p);

     }
     runTime.write();
     Info << "ExecutionTime = " << runTime.elapsedCpuTime()   
          << " s" << "ClockTime = "<< runTime.elapsedClockTime() 
          << " s" << nl << endl;

  }

    Info<< "\nEnd\n" << endl;
    return 0;
}

// ************************************************************************* //
