/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    engineSwirl

Description
    Generates a swirling flow for engine calulations.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    forAll(mesh.C(), celli)
    {
        vector location = mesh.C()[celli];
        if (location[1] > 0.45){
            c[celli] = 1;
        } 
        else if (location[1] > 0.43 && location[1] <= 0.45){
            c[celli] = 0.9;
        }
        else {
            c[celli] = 0;
        }
        //Info << location << endl;
    }

    //c.correctBoundaryConditions();
    c.write();

    Info<< "\n end\n";

    return 0;
}


// ************************************************************************* //
