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

        volVectorField dpCont = 2.0 / 3.0 * fvc::grad(phase1.turbulence().k() * alpha1);
        if(!KM.dp())
        {
            dpCont *= 0.0;
        }

        U1Eqn =
        (
            fvm::ddt(alpha1, U1)
          + fvm::div(alphaPhi1, U1)
          - fvm::Sp(fvc::ddt(alpha1) + fvc::div(alphaPhi1), U1)
          + phase1.turbulence().divDevReff(U1)
          + dpCont
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
        if(pimple.momentumPredictor())
        {
            solve(U1Eqn == - alpha1 * fvc::grad(p) / rho1);
            //fvOptions.correct(U1);
        }
    }

    {

        //volScalarField k = fluid.otherPhase(KM.dispersedPhase()).turbulence().k();
        //volScalarField epsilon = fluid.otherPhase(KM.dispersedPhase()).turbulence().epsilon();
        volScalarField epsilon = KM.dispersedPhase().turbulence().epsilon();
        volScalarField k = KM.dispersedPhase().turbulence().k();
        //dimensionedScalar maxT("maxT", k.dimensions(), 5);
        //dimensionedScalar maxEps("maxT", epsilon.dimensions(), 10);
        //k = min(k, maxT);
        //epsilon = min(epsilon, maxEps);

        //volScalarField k =
            //min(
                //k0,
                //pow(mag(KM.dispersedPhase().U()), 2)
            //);
        //k.boundaryField() = k0.boundaryField();

	//-------------------------turbulent only
	int updateCount = 5;
	for(int i = 0; i < updateCount; i++) KM.update(k, epsilon, 0);

        Info << "Calculating pressure correction" << endl;
	volScalarField X2 = KM.pressureCorrection();
	//volVectorField dp(2.0 / 3.0 * fvc::grad(X2 * k * alpha2));
        if(KM.dp())
        {
            dp = 2.0 / 3.0 * fvc::grad(X2 * k * alpha2);
            dp.correctBoundaryConditions();
        }else
        {
            dp = dp * 0.0;
        }
	
	volVectorField& F(KM.collisionalF(alphaPhi2));
	volScalarField FSp(KM.collisionalSp(alphaPhi2));
        KM.print();

        F_old = F;

        U2Eqn =
        (
            fvm::ddt(alpha2, U2)
          + fvm::div(alphaPhi2, U2)
          - fvm::Sp(fvc::ddt(alpha2) + fvc::div(alphaPhi2), U2)
          //+ phase2.turbulence().divDevReff(U2)
          + KM.divDevReff(U2)
         ==
          fvm::Sp(FSp, U2)
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

        if(pimple.momentumPredictor())
        {
            solve(U2Eqn == - alpha2 * fvc::grad(p) / rho2);
            //fvOptions.correct(U2);
        }
    }


}



