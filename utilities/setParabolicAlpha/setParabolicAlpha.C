/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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
    twoPhaseEulerFoam

Description
    Solver for a system of 2 compressible fluid phases with one phase
    dispersed, e.g. gas bubbles in a liquid including heat-transfer.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "twoPhaseSystem.H"
#include "PhaseIncompressibleTurbulenceModel.H"
#include "pimpleControl.H"
#include "IOMRFZoneList.H"
#include "fixedFluxPressureFvPatchScalarField.H"
#include "kineticFluidModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noParallel();
    timeSelector::addOptions();
    #include "setRootCase.H"

    #include "createTime.H"
    #include "createMesh.H"
    #include "readGravitationalAcceleration.H"
    #include "createFields.H"
    #include "createMRFZones.H"
    #include "initContinuityErrs.H"
    #include "readTimeControls.H"
    #include "CourantNos.H"
    #include "setInitialDeltaT.H"

    pimpleControl pimple(mesh);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    instantList timeDirs = timeSelector::select0(runTime, args);
    runTime.setTime(timeDirs[0], 0);

    Info<< "\nStarting time loop\n" << endl;

    forAll(timeDirs, timeI)
    {

        scalar h = 0.076;
        scalar x0 = h / 2;
        scalar alpha0 = 0.4;
        scalar alpha_h = 0.25;

        dimensionedScalar dL("dL", dimLength, 1);
        volScalarField x = mesh.C().component(0) / dL;
        //volScalarField y = mesh.C().component(1);

        volScalarField r = x - x0;//sqrt(pow(x, 2) + pow(y, 2));
        
        dimensionedScalar a = 3.0 * (alpha0 - alpha_h) / pow(h / 2, 2);
        dimensionedScalar b = 3.0 * alpha0 - 2.0 * alpha_h;

        volScalarField A
            (
                "A",
                - a * pow(r, 2) + b
            );
        volScalarField B("B", 1.0 - A);
        A.write();
        B.write();

    }


    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
