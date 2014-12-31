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
          //+ 2.0 / 3.0 * fvc::grad(phase1.turbulence().k() * alpha1)
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

        //volScalarField k = fluid.otherPhase(KM.dispersedPhase()).turbulence().k();
        //volScalarField epsilon = fluid.otherPhase(KM.dispersedPhase()).turbulence().epsilon();
        volScalarField epsilon = KM.dispersedPhase().turbulence().epsilon();
        volScalarField k = KM.dispersedPhase().turbulence().k();
        dimensionedScalar maxT("maxT", k.dimensions(), 10);
        dimensionedScalar maxEps("maxT", epsilon.dimensions(), 100);
        k = min(k, maxT);
        epsilon = min(epsilon, maxEps);

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
	volVectorField dp(2.0 / 3.0 * fvc::grad(X2 * k * alpha2));
	
        Info << "Calculating F1" << endl;
	volVectorField F1 ("F1", KM.F1(alphaPhi2));
        Info << "Calculating F2" << endl;
	volVectorField F2 ("F2", KM.F2(alphaPhi2));
        Info << "Calculating F3" << endl;
	volVectorField F3 ("F3", KM.F3(alphaPhi2));


	F = F1 + F3;
	volVectorField F0 ("0", F);
        dimensionedScalar F_tres ("Ftres", dimLength / pow(dimTime, 2), 1000.0);
        //F.boundaryField() = vector(0, 0, 0);
        //F0.boundaryField() = vector(0, 0, 0);
        //F1.boundaryField() = vector(0, 0, 0);
        //F3.boundaryField() = vector(0, 0, 0);
        //F2.boundaryField() = vector(0, 0, 0);
        /*F.component(0) = min(F.component(0), F_tres);;*/
        //F.component(0) = max(F.component(0), -F_tres);;
        //F.component(1) = min(F.component(1), F_tres);;
        //F.component(1) = min(F.component(1), F_tres);;
        //F.component(2) = max(F.component(2), -F_tres);;
        //F.component(2) = max(F.component(2), -F_tres);;
        
        //const fvPatchList& patches = mesh.boundary();

        //forAll(patches, patchi)
        //{
                //const fvPatch& curPatch = patches[patchi];
                //if (!isType<wallFvPatch>(curPatch))
                //{
                        //forAll(curPatch, facei)
                        //{
                                //label faceCelli = curPatch.faceCells()[facei];
                                //F0[faceCelli] = vector(0, 0, 0);
                                //F[faceCelli] = vector(0, 0, 0);
                        //}
                //}
        //}
        //forAll(mesh.C(), celli)
        //{
                //F0[celli].x() = min(F0[celli].x(), 500.0);
                //F0[celli].y() = min(F0[celli].y(), 500.0);
                //F0[celli].z() = min(F0[celli].z(), 500.0);

                //F0[celli].x() = max(F0[celli].x(), -500.0);
                //F0[celli].y() = max(F0[celli].y(), -500.0);
                //F0[celli].z() = max(F0[celli].z(), -500.0);

                //F[celli].x() = min(F[celli].x(), 500.0);
                //F[celli].y() = min(F[celli].y(), 500.0);
                //F[celli].z() = min(F[celli].z(), 500.0);

                //F[celli].x() = max(F[celli].x(), -500.0);
                //F[celli].y() = max(F[celli].y(), -500.0);
                //F[celli].z() = max(F[celli].z(), -500.0);
        //}
        forAll(mesh.C(), celli)
        {
                //scalar count = 1.0;
                //forAll(mesh.cellCells()[celli], cellj)
                //{
                    //label n_id = mesh.cellCells()[celli][cellj];
                    //count += 1.0;
                    //F[celli] += F0[n_id];
                //}
                //F[celli] /= count;

                //if(mesh.C()[celli].z() < 0.2 || mesh.C()[celli].z() > 3.0)
                //{
                    //F[celli] = vector(0, 0, 0);
                //}
                //F[celli].z() = 0.0f;
                //F[celli].x() = min(F[celli].x(), 1000.0);
                //F[celli].y() = min(F[celli].y(), 1000.0);
                //F[celli].z() = min(F[celli].z(), 1000.0);

                //F[celli].x() = max(F[celli].x(), -1000.0);
                //F[celli].y() = max(F[celli].y(), -1000.0);
                //F[celli].z() = max(F[celli].z(), -1000.0);
        }
        //F = 0.5 * (F + F_old);



	
	// deleta force next to wall
        //const fvPatchList& patches = mesh.boundary();

        //forAll(patches, patchi)
        //{
                //const fvPatch& curPatch = patches[patchi];
                //if (isType<wallFvPatch>(curPatch))
                //{
                        //forAll(curPatch, facei)
                        //{
                                //label faceCelli = curPatch.faceCells()[facei];
                                //F[faceCelli] = vector(0, 0, 0);
                        //}
                //}
        //}
	//-------------------------

        F.correctBoundaryConditions();
	volScalarField Fmag = mag(F);
        if(mesh.time().outputTime())
        {
            F1.write();
            F2.write();
            F3.write();
            //F.write();
        }
        F_old = F;

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
	Info << "Fmag: " << Fmag.weightedAverage(mesh.V()).value()
	  <<" min: " << min(Fmag).value()
	  <<" max: " << max(Fmag).value()
	  <<" dims: " << Fmag.dimensions()
	  << endl;
	Info << "------------------------------------" << endl;

        U2Eqn =
        (
            fvm::ddt(alpha2, U2)
          + fvm::div(alphaPhi2, U2)
          - fvm::Sp(fvc::ddt(alpha2) + fvc::div(alphaPhi2), U2)
          + phase2.turbulence().divDevReff(U2)
          //+ KM.divDevReff(U2)
         ==
          F 
	  //- dp 
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
