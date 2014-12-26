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
        //volScalarField epsilon = KM.dispersedPhase().turbulence().epsilon();
        //volScalarField k = KM.dispersedPhase().turbulence().k();
        volScalarField k = fluid.otherPhase(KM.dispersedPhase()).turbulence().k();
        volScalarField epsilon = fluid.otherPhase(KM.dispersedPhase()).turbulence().epsilon();

        int updateCount = 5;
        for(int i = 0; i < updateCount; i++) KM.update(k, epsilon, 0);

        volScalarField cd
            (
                "cd",
                fluid.dragCoeff() * KM.dispersedPhase() / KM.dispersedPhase().rho()
            );


        volScalarField pressureCorrFactor = KM.pressureCorrection();

        volVectorField dp
            (
                "dp",
                2.0 / 3.0 
                * fvc::grad
                (
                    phase2.turbulence().k() 
                    * pressureCorrFactor
                    * KM.dispersedPhase()
                )
            );


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
                "tau",
                volScalarField(KM.tauTotal())
            );

        cd.write();
        dp.write();
        ratio.write();
        correctedViscosity.write();
        volScalarField X("X", pressureCorrFactor);
        X.write();
        tau.write();

        volScalarField J1("J1", KM.J1());
        volScalarField J2("J2", KM.J2());
        volScalarField J3("J3", KM.J3());
        volScalarField J4("J4", KM.J4());
        volScalarField beta1("beta_1", KM.beta1());
        volScalarField beta2("beta_2", KM.beta2());
        J1.write();
        J2.write();
        J3.write();
        J4.write();
        beta1.write();
        beta2.write();
    }


    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
