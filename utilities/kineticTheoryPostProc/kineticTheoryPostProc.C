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
      //--------------turblent corrections
      runTime.setTime(timeDirs[timeI], timeI);
      volScalarField epsilon = KM.dispersedPhase().turbulence().epsilon();
      volScalarField k = KM.dispersedPhase().turbulence().k();
      KM.update(k, epsilon);

      volScalarField cd
      (
          "cd",
          fluid.dragCoeff() * KM.dispersedPhase() / KM.dispersedPhase().rho()
      );


      volScalarField pressureCorrFactor = KM.pressureCorrection();

      volVectorField dp
        (
            "dp_t",
            2.0 / 3.0 
            * fvc::grad
            (
                phase2.turbulence().k() 
                * pressureCorrFactor
                * KM.dispersedPhase()
            )
        );


    //volScalarField k =  KM.dispersedPhase().turbulence().k();
    volScalarField correctedViscosity
        (
            "dNut",
            8.0 / 9.0 
            * KM.a() * (1.0 / (KM.E1() + 2.0 * k * KM.E2())) * pow(k, 2)
        );
    
    volScalarField ratio
        (
            "nutCorrRatio",
            0.5 * correctedViscosity 
            / KM.dispersedPhase().turbulence().nuEff()
        );

    volScalarField tau
        (
            "tau_t",
            volScalarField(KM.tau())
        );

      cd.write();
      dp.write();
      ratio.write();
      correctedViscosity.write();
      volScalarField X("X_t", pressureCorrFactor);
      X.write();
      tau.write();
      //--------------laminar corrections
      volScalarField p_mod("p_mod", p / rho2 - g.component(2) * mesh.C().component(2));
      volScalarField Umag2 = 0.5 * pow(mag(U2), 2) 
          + dimensionedScalar("minU", pow(dimLength,2) / pow(dimTime,2), 0.8);
      volScalarField dissipation = Umag2 * cd;
      KM.update(Umag2, dissipation);

      volScalarField pressureCorrFactor2 = KM.pressureCorrection();

      volVectorField dp2
        (
            "dp_l",
            2.0 / 3.0 
            * fvc::grad
            (
                p_mod
                * pressureCorrFactor2
                * KM.dispersedPhase()
            )
        );

    volScalarField tau2
        (
            "tau_l",
            volScalarField(KM.tau())
        );

      dp2.write();
      volScalarField X2("X_l", pressureCorrFactor2);
      X2.write();
      tau2.write();
      p_mod.write();
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
