mrfZones.correctBoundaryVelocity(U1);
mrfZones.correctBoundaryVelocity(U2);
mrfZones.correctBoundaryVelocity(U);

fvVectorMatrix U1Eqn(U1, U1.dimensions()*dimVol/dimTime);
fvVectorMatrix U2Eqn(U2, U2.dimensions()*dimVol/dimTime);

volScalarField dragCoeff(fluid.dragCoeff());

{
    volScalarField virtualMassCoeff(fluid.virtualMassCoeff());

    volVectorField liftForce(fluid.liftForce());
    volVectorField wallLubricationForce(fluid.wallLubricationForce());
    volVectorField turbulentDispersionForce(fluid.turbulentDispersionForce());

    volScalarField pressureCorrFactor = KM.pressureCorrection();
    Info << "pressure correction: " << pressureCorrFactor.weightedAverage(mesh.V()).value()
      <<" min: " << min(pressureCorrFactor).value()
      <<" max: " << max(pressureCorrFactor).value() << endl;

    {
	
	volScalarField X = pressureCorrFactor;
	if(&alpha1 == &KM.dispersedPhase())
	{
	  Info << "correcting pressure in phase1 " << endl;
	  Info << "phase name is" << phase1.name()
	    <<" and dispersed phase name is " << KM.dispersedPhase().name() << endl;

	  volVectorField dp = fvc::grad(p * X) / rho1;
	  Info << "pressure correction: " << dp.weightedAverage(mesh.V()).value()
	    <<" min: " << min(dp).value()
	    <<" max: " << max(dp).value() << endl;

	}else
	{
	  X = 0.0 * X;
	}

        U1Eqn =
        (
            fvm::ddt(alpha1, U1)
          + fvm::div(alphaPhi1, U1)
          - fvm::Sp(fvc::ddt(alpha1) + fvc::div(alphaPhi1), U1)
          + phase1.turbulence().divDevReff(U1)
         ==
          - fvm::Sp(dragCoeff/rho1, U1)
          - alpha1*alpha2/rho1
           *(
                liftForce
              + wallLubricationForce
              + turbulentDispersionForce
            )
          - virtualMassCoeff/rho1
           *(
                fvm::ddt(U1)
              + fvm::div(phi1, U1)
              - fvm::Sp(fvc::div(phi1), U1)
              - DDtU2
            )
        );
        mrfZones.addCoriolis(alpha1 + virtualMassCoeff/rho1, U1Eqn);
        U1Eqn.relax();
    }

    {

	volScalarField X = pressureCorrFactor;
	if(&alpha2 == &KM.dispersedPhase())
	{
	  Info << "correcting pressure in phase 2" << endl;
	  Info << "phase name is" << phase2.name()
	    <<" and dispersed phase name is " << KM.dispersedPhase().name() << endl;
	}else
	{
	  X = 0.0 * X;
	}

        //volVectorField dp = alpha2 * 2.0 / 3.0  * fvc::grad(phase2.turbulence().k() * X);
	//volVectorField dp = alpha2 * X fvc::grad(X * 2.0 / 3.0 * phase2.turbulence().k());
        //volVectorField dp = X * fvc::grad(p) / rho2 
            //+ 2.0 / 3.0 * phase2.turbulence().k() * fvc::grad(X);
        volVectorField dp = 2.0 / 3.0 * fvc::grad(phase2.turbulence().k() * X);
	Info << "pressure correction: " << dp.weightedAverage(mesh.V()).value()
	  <<" min: " << min(dp).value()
	  <<" max: " << max(dp).value() << endl;

        U2Eqn =
        (
            fvm::ddt(alpha2, U2)
          + fvm::div(alphaPhi2, U2)
	  + dp * alpha2
          - fvm::Sp(fvc::ddt(alpha2) + fvc::div(alphaPhi2), U2)
          //+ phase2.turbulence().divDevReff(U2)
          + KM.divDevReff(U2)
         ==
          - fvm::Sp(dragCoeff/rho2, U2)
          + alpha1*alpha2/rho2
           *(
                liftForce
              + wallLubricationForce
              + turbulentDispersionForce
            )
          - virtualMassCoeff/rho2
           *(
                fvm::ddt(U2)
              + fvm::div(phi2, U2)
              - fvm::Sp(fvc::div(phi2), U2)
              - DDtU1
            )
        );
        mrfZones.addCoriolis(alpha2 + virtualMassCoeff/rho2, U2Eqn);
        U2Eqn.relax();
    }
}
