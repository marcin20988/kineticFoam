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
	/*volScalarField X = pressureCorrFactor;
	if(&alpha2 == &KM.dispersedPhase())
	{
	  Info << "correcting pressure in phase 2" << endl;
	  Info << "phase name is" << phase2.name()
	    <<" and dispersed phase name is " << KM.dispersedPhase().name() << endl;
	}else
	{
	  X = 0.0 * X;
	}*/
	volScalarField p_mod = p / rho2 - g.component(2) * mesh.C().component(2);
        Info << "Updating laminar corrections " << endl;
        //KM.update(3.0 * p_mod / 2.0);
	volScalarField Umag2 = 0.5 * pow(mag(U2), 2) 
		+ dimensionedScalar("minU", pow(dimLength,2) / pow(dimTime,2), 0.2);
        KM.update(Umag2);
        //KM.update(3.0 * p_mod / 2.0);
        Info << "-------------------------------" << endl;
	volScalarField X2 = KM.pressureCorrection();
        //KM.update(p_mod);
	volVectorField dp2 = fvc::grad(p_mod * (X2 - 1.0) * alpha2) ;

        Info << "Updating turbulent corrections " << endl;
        KM.update(k);
        Info << "-------------------------------" << endl;
	volScalarField X = KM.pressureCorrection();
        volVectorField dp = 2.0 / 3.0 * fvc::grad(phase2.turbulence().k() * alpha2);


        //KM.update(KM.dispersedPhase().turbulence().k());

	Info << "modified pressure: " << p_mod.weightedAverage(mesh.V()).value()
	  <<" min: " << min(p_mod).value()
	  <<" max: " << max(p_mod).value() << endl;
	Info << "turbulent correction: " << dp.weightedAverage(mesh.V()).value()
	  <<" min: " << min(dp).value()
	  <<" max: " << max(dp).value() << endl;
	Info << "kinetic pressure correction: " << dp2.weightedAverage(mesh.V()).value()
	  <<" min: " << min(dp2).value()
	  <<" max: " << max(dp2).value() << endl;
	Info << "Umag2: " << Umag2.weightedAverage(mesh.V()).value()
	  <<" min: " << min(Umag2).value()
	  <<" max: " << max(Umag2).value() << endl;

        U2Eqn =
        (
            fvm::ddt(alpha2, U2)
          + fvm::div(alphaPhi2, U2)
          - fvm::Sp(fvc::ddt(alpha2) + fvc::div(alphaPhi2), U2)
          + phase2.turbulence().divDevReff(U2)
          //+ KM.divDevReff(U2)
         ==
	  - dp 
	  - dp2
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
