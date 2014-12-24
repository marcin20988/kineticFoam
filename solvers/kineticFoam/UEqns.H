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

    /*KM.update(k);
    volScalarField pressureCorrFactor = KM.pressureCorrection();
    Info << "pressure correction: " << pressureCorrFactor.weightedAverage(mesh.V()).value()
      <<" min: " << min(pressureCorrFactor).value()
      <<" max: " << max(pressureCorrFactor).value() << endl;*/

    {
	
	/*volScalarField X = pressureCorrFactor;
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
	}*/

        U1Eqn =
        (
            fvm::ddt(alpha1, U1)
          + fvm::div(alphaPhi1, U1)
          - fvm::Sp(fvc::ddt(alpha1) + fvc::div(alphaPhi1), U1)
          + phase1.turbulence().divDevReff(U1)
          + 2.0 / 3.0 * fvc::grad(phase1.turbulence().k() * alpha1)
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

        volScalarField k = KM.dispersedPhase().turbulence().k();
        volScalarField epsilon = KM.dispersedPhase().turbulence().epsilon();

	//-------------------------turbulent only
	int updateCount = 10;
	for(int i = 0; i < updateCount; i++) KM.update(k, epsilon, 0);
	volScalarField X2 = KM.pressureCorrection();
	volVectorField dp(2.0 / 3.0 * fvc::grad(X2 * k * alpha2));
	
	volVectorField F1 = KM.F1(alphaPhi2);
	volVectorField F2 = KM.F2(alphaPhi2);
	volVectorField F3 = KM.F3(alphaPhi2);
	volVectorField F = F1 + F2 + F3;

	Info << "------------------------------------" << endl;
	Info << "correction factor: " << X2.weightedAverage(mesh.V()).value()
	  <<" min: " << min(X2).value()
	  <<" max: " << max(X2).value() << endl;
	Info << "force: " << dp.weightedAverage(mesh.V()).value()
	  <<" min: " << min(dp).value()
	  <<" max: " << max(dp).value() << endl;
	Info << "F1: " << F1.weightedAverage(mesh.V()).value()
	  <<" min: " << min(F1).value()
	  <<" max: " << max(F1).value() 
	  <<" dims: " << F1.dimensions()
	  << endl;
	Info << "F2: " << F2.weightedAverage(mesh.V()).value()
	  <<" min: " << min(F2).value()
	  <<" max: " << max(F2).value()
	  <<" dims: " << F2.dimensions()
	  << endl;
	Info << "F3: " << F3.weightedAverage(mesh.V()).value()
	  <<" min: " << min(F3).value()
	  <<" max: " << max(F3).value()
	  <<" dims: " << F3.dimensions()
	  << endl;
	Info << "F: " << F.weightedAverage(mesh.V()).value()
	  <<" min: " << min(F).value()
	  <<" max: " << max(F).value()
	  <<" dims: " << F.dimensions()
	  << endl;
	Info << "------------------------------------" << endl;

        U2Eqn =
        (
            fvm::ddt(alpha2, U2)
          + fvm::div(alphaPhi2, U2)
          - fvm::Sp(fvc::ddt(alpha2) + fvc::div(alphaPhi2), U2)
          //+ phase2.turbulence().divDevReff(U2)
          + KM.divDevReff(U2)
         ==
	  - dp 
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
