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

\*---------------------------------------------------------------------------*/

#include "kineticFluidModel.H"


namespace Foam
{
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

kineticFluidModel::kineticFluidModel(const twoPhaseSystem& fluid):
  fluid_(fluid),
  mesh_(fluid_.mesh()),
  kineticFluidModelDict_
  (
   IOdictionary
   (
    IOobject
    (
     "kineticFluidModelProperties",
     fluid_.mesh().time().constant(),
     fluid_.mesh(),
     IOobject::MUST_READ_IF_MODIFIED,
     IOobject::NO_WRITE
    )
   )
  ),
  continuousPhaseName_(kineticFluidModelDict_.lookup("continuousPhase")),
  dispersedPhaseName_(kineticFluidModelDict_.lookup("dispersedPhase")),
  tau_
  (
   relaxationTime::New
   (
     fluid_,
     kineticFluidModelDict_,
     dispersedPhaseName_,
     dispersedPhase(),
     *this
   )
  ),
  E1_
  (
   IOobject
   (
    "E1",
    mesh_.time().timeName(),
    mesh_,
    IOobject::MUST_READ,
    IOobject::AUTO_WRITE
   ),
   mesh_
   //dimensionedScalar("E1", dimless, 0.0)
  ),
  E2_
  (
   IOobject
   (
    "E2",
    mesh_.time().timeName(),
    mesh_,
    IOobject::MUST_READ,
    IOobject::AUTO_WRITE
   ),
   mesh_
   //dimensionedScalar("E2", pow(dimLength, -2) / pow(dimTime, -2), 0.0)
  ),
  a_
  (
   IOobject
   (
    "E2",
    mesh_.time().timeName(),
    mesh_,
    IOobject::NO_READ,
    IOobject::AUTO_WRITE
   ),
   mesh_,
   dimensionedScalar("a", pow(dimLength, -2) / pow(dimTime, -3), 0.0)
  ),
  T_
  (
   IOobject
   (
    "kineticTemperature",
    mesh_.time().timeName(),
    mesh_,
    IOobject::NO_READ,
    IOobject::AUTO_WRITE
   ),
   mesh_,
   dimensionedScalar("a", pow(dimLength, 2) / pow(dimTime, 2), 0.0)
  ),
  epsilon_
  (
   IOobject
   (
    "kineticTDissipation",
    mesh_.time().timeName(),
    mesh_,
    IOobject::NO_READ,
    IOobject::AUTO_WRITE
   ),
   mesh_,
   dimensionedScalar("a", pow(dimLength, 2) / pow(dimTime, 3), 0.0)
  ),
  deltaG_
  (
   IOobject
   (
    "deltaG_",
    mesh_.time().timeName(),
    mesh_,
    IOobject::NO_READ,
    IOobject::AUTO_WRITE
   ),
   mesh_,
   dimensionedVector("deltaG", dimless / dimLength, vector(0.0, 0.0, 0.0))
  ),
  e_(readScalar(kineticFluidModelDict_.lookup("e"))),
  maxF_(readScalar(kineticFluidModelDict_.lookup("maxF")))
{}



// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<kineticFluidModel>
kineticFluidModel::New(const twoPhaseSystem& fluid)
{
    return autoPtr<kineticFluidModel>(new kineticFluidModel(fluid));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const phaseModel& kineticFluidModel::dispersedPhase() const
{
  if(fluid_.phase1().name() == dispersedPhaseName_)
  {
    return fluid_.phase1();
  }
  else if(fluid_.phase2().name() == dispersedPhaseName_)
  {
    return fluid_.phase2();
  }
  else
  {
    FatalErrorIn("phaseModel::dispersedPhae")
      << "there is no phase: " << dispersedPhaseName_ << endl
      << exit(FatalError);
    return fluid_.phase2();
  }
};

// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

void kineticFluidModel::operator=(const kineticFluidModel& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorIn("Foam::kineticFluidModel::operator=(const Foam::kineticFluidModel&)")
            << "Attempted assignment to self"
            << abort(FatalError);
    }
}

void kineticFluidModel::update
(
    const volScalarField& T,
    const volScalarField& epsilon,
    const int check
)
{
  Info << "Updating kinetic model" << endl;
  T_ = T;
  epsilon_ = epsilon;
  volScalarField k =  T;//dispersedPhase().turbulence().k();
  //volScalarField epsilon =  dispersedPhase().turbulence().epsilon();

  volScalarField tau = tau_ -> total();
  if(check == 0)
  {
	tau = tau_ -> turbulent();
  }else if(check == 1)
  {
	tau = tau_ -> laminar();
  }
  // in the derivation drag coefficient is taken in units of s^-1
  // and is strictly negative
  volScalarField cd = - fluid_.dragCoeff() 
      / dispersedPhase().rho() * dispersedPhase();

  volScalarField cd_stokes = - 3.0 * 3.14 * fluid_.otherPhase(dispersedPhase()).nu()
      * fluid_.otherPhase(dispersedPhase()).d() 
      / (3.13 / 6.0 * pow(fluid_.otherPhase(dispersedPhase()).d(), 3));





  Info << "Collisional relaxation time: " 
    << tau.weightedAverage(tau.mesh().V()).value()
    <<" min: " << min(tau).value()
    <<" max: " << max(tau).value() << endl;
  Info << "Stokes drag: " << cd_stokes.weightedAverage(tau.mesh().V()).value() << endl;
  Info << "Drag coefficient: " << cd.weightedAverage(tau.mesh().V()).value()
    <<" min: " << min(cd).value()
    <<" max: " << max(cd).value() << endl;
  //Info << "using stokes drag" << endl;
  //cd = cd_stokes;

  E1_ = 
    (
        2.0 * k - 8 * tau * cd * k 
        + 3.0 * tau * epsilon - 6.0 * cd * pow(tau, 2) * epsilon
    )
    /
    (
        16.0 * pow(cd, 2) * pow(tau, 2) * k - 12.0 * cd * tau * k + 2.0 * k
    );


  //NOTE: this is E2 * k !!!
  E2_ = 3.0 * tau * epsilon / (4.0 * k * (1.0 - 4.0 * tau * cd));
  E2_ = min(1.0, E2_);

  E1_.correctBoundaryConditions();
  E2_.correctBoundaryConditions();
  //E2_.boundaryField() = 0;


  Info << "E1: " << E1_.weightedAverage(tau.mesh().V()).value()
    <<" min: " << min(E1_).value()
    <<" max: " << max(E1_).value() << endl;
  Info << "E2: " << E2_.weightedAverage(tau.mesh().V()).value()
    <<" min: " << min(E2_).value()
    <<" max: " << max(E2_).value() << endl;
  
  a_ = 2.0 * tau / (k *(1.0 - 6.0 * tau * cd));
  Info << "a: " << a_.weightedAverage(tau.mesh().V()).value()
    <<" min: " << min(a_).value()
    <<" max: " << max(a_).value() << endl;

  deltaG_ = deltaG();
  //deltaG_.boundaryField() = vector(0, 0, 0);

  scalar treshold = 1e04;
  //deltaG_.component(0) =
      //min(treshold, deltaG_.component(0));
  //deltaG_.component(0) =
      //max(-treshold, deltaG_.component(0));
  forAll(deltaG_.mesh().V(), celli)
  {
      deltaG_[celli].x() = min(treshold, deltaG_[celli].x());
      deltaG_[celli].x() = max(-treshold, deltaG_[celli].x());

      deltaG_[celli].y() = min(treshold, deltaG_[celli].y());
      deltaG_[celli].y() = max(-treshold, deltaG_[celli].y());

      deltaG_[celli].z() = min(treshold, deltaG_[celli].z());
      deltaG_[celli].z() = max(-treshold, deltaG_[celli].z());
  }

  //deltaG_.internalField().component(1) =
      //min(1e04, deltaG_.internalField().component(1));
  //deltaG_.internalField().component(1) =
      //max(-1e04, deltaG_.internalField().component(1));

  //deltaG_.internalField().component(2) =
      //min(1e04, deltaG_.internalField().component(2));
  //deltaG_.internalField().component(2) =
      //max(-1e04, deltaG_.internalField().component(2));

  //const fvPatchList& patches = deltaG_.mesh().boundary();

  //forAll(patches, patchi)
  //{
      //const fvPatch& curPatch = patches[patchi];
      //if (isType<wallFvPatch>(curPatch))
      //{
          //forAll(curPatch, facei)
          //{
              //label faceCelli = curPatch.faceCells()[facei];
              //deltaG_[faceCelli] = vector(0, 0, 0);
          //}
      //}
  //}

  Info << "deltaG: " << deltaG_.weightedAverage(tau.mesh().V()).value()
      <<" min: " << min(deltaG_).value()
      <<" max: " << max(deltaG_).value() << endl;
};


tmp<volScalarField> kineticFluidModel::pressureCorrection() const
{
  volScalarField k =  temp();//dispersedPhase().turbulence().k();

  return (1.0 + E1_ + 10.0 / 3.0 * E2_) / (1.0 + E1_ + 2.0 * E2_);
};


tmp<fvVectorMatrix> kineticFluidModel::divDevReff(const volVectorField& U)
{
    volScalarField k =  temp();

    volScalarField correctedViscosity = 0.5 * 8.0 / 9.0 
        * a_ * (1.0 / (1.0 + E1_ + 2.0 * E2_)) * pow(k, 2);
    
    correctedViscosity += dispersedPhase().nu();

    volScalarField ratio = correctedViscosity 
        / dispersedPhase().turbulence().nuEff();
    Info << "Ratio of viscosity correction: " << ratio.weightedAverage(k.mesh().V()).value()
        <<" min: " << min(ratio).value()
        <<" max: " << max(ratio).value() << endl;

    return
        (
            - fvm::laplacian(correctedViscosity, U)
            - fvc::div(correctedViscosity*dev(T(fvc::grad(U))))
        );
}


tmp<volScalarField> kineticFluidModel::tauTurbulent()
{
    return tau_ -> turbulent();
}

tmp<volScalarField> kineticFluidModel::tauLaminar()
{
    return tau_ -> laminar();
}

tmp<volScalarField> kineticFluidModel::tauTotal()
{
    return tau_ -> total();
}

tmp<volScalarField> kineticFluidModel::beta(scalar c1, scalar c2, scalar c3) const
{
	volScalarField F1 = 1 + E1_;
	return c1 * pow(F1, 2) + c2 * F1 * E2_ + c3 * pow(E2_,2);
}
tmp<volScalarField> kineticFluidModel::beta1() const
{
        return beta(0.14, 1.51, 2.15);
}

tmp<volScalarField> kineticFluidModel::beta2() const
{
        return beta(0.23, 1.81, 2.63);
}

tmp<volScalarField> kineticFluidModel::beta3() const
{
        return beta(0.21, 2.26, 3.22);
}

tmp<volScalarField> kineticFluidModel::beta4() const
{
        return beta(0.57, 4.53, 6.57);
}

tmp<volScalarField> kineticFluidModel::beta5() const
{
        return beta(-0.10, -1.68, -2.61);
}

tmp<volScalarField> kineticFluidModel::beta6() const
{
        return beta(-0.01, -0.32, -0.57);
}

tmp<volScalarField> kineticFluidModel::beta7() const
{
        return beta(0.08, 0.64, 0.87);
}

tmp<volScalarField> kineticFluidModel::beta8() const
{
        return beta(0.09, 1.49, 2.32);
}

tmp<volScalarField> kineticFluidModel::beta9() const
{
        return beta(-0.05, -0.71, -1.04);
}

tmp<volScalarField> kineticFluidModel::beta10() const
{
        return beta(0.05, 0.57, 0.80);
}

tmp<volScalarField> kineticFluidModel::beta11() const
{
        return beta(0.14, 1.13, 1.64);
}

tmp<volScalarField> kineticFluidModel::beta12() const
{
        return beta(0.06, 0.80, 1.17);
}

tmp<volScalarField> kineticFluidModel::beta13() const
{
        return beta(0.24, 2.23, 3.21);
}

tmp<volScalarField> kineticFluidModel::beta14() const
{
        return beta(0.0, 0.0, 0.53);
}

tmp<volScalarField> kineticFluidModel::beta15() const
{
        return beta(-0.05, -0.71, -1.04);
}

tmp<volScalarField> kineticFluidModel::beta16() const
{
        return beta(0.11, 1.13, 1.61);
}

tmp<volScalarField> kineticFluidModel::beta17() const
{
        return beta(0.42, 3.40, 4.93);
}

tmp<volScalarField> kineticFluidModel::J1() const
{
	dimensionedScalar smallT("smallT", T_.dimensions(), 1e-08);
	volScalarField a = min
            ( 
                mag(dispersedPhase().U()) * sqrt(3.0 / (8.0 * T_ + smallT)),
                20.0
            );
        a = max(a, 0.1);

	Info << "a-function coefficient: " 
		<< a.weightedAverage(a.mesh().V()).value()
		<<" min: " << min(a).value()
		<<" max: " << max(a).value() << endl;
        volScalarField j1 = beta1() + pow(a, 2) * beta2();
        //j1.boundaryField() = 0;
        // deleta force next to wall
        //const fvPatchList& patches = j1.mesh().boundary();

        //forAll(patches, patchi)
        //{
            //const fvPatch& curPatch = patches[patchi];
            //if (isType<wallFvPatch>(curPatch))
            //{
                //forAll(curPatch, facei)
                //{
                    //label faceCelli = curPatch.faceCells()[facei];
                    //j1[faceCelli] = 0;
                //}
            //}
        //}
	//-------------------------
	return - T_ * j1;
}

tmp<volScalarField> kineticFluidModel::J2() const
{
	dimensionedScalar smallT("smallT", T_.dimensions(), 1e-08);
	//volScalarField a = mag(dispersedPhase().U()) * sqrt(3.0 / (8.0 * T_));
	volScalarField a = min
            ( 
                mag(dispersedPhase().U()) * sqrt(3.0 / (8.0 * T_ + smallT)),
                20.0
            );
        a = max(a, 0.1);

        volScalarField j2 = beta3() + pow(a, 2) * beta4();
        //j2.boundaryField() = 0;
        // deleta force next to wall
        const fvPatchList& patches = j2.mesh().boundary();

        //forAll(patches, patchi)
        //{
            //const fvPatch& curPatch = patches[patchi];
            //if (isType<wallFvPatch>(curPatch))
            //{
                //forAll(curPatch, facei)
                //{
                    //label faceCelli = curPatch.faceCells()[facei];
                    //j2[faceCelli] = 0;
                //}
            //}
        //}
	//-------------------------
	return - T_ * j2;
}

tmp<volScalarField> kineticFluidModel::J3() const
{
	dimensionedScalar smallT("smallT", T_.dimensions(), 1e-08);
	//volScalarField a = max( mag(dispersedPhase().U()) * sqrt(3.0 / (8.0 * T_ + smallT)), 0.5);
	volScalarField a = min
            ( 
                mag(dispersedPhase().U()) * sqrt(3.0 / (8.0 * T_ + smallT)),
                20.0
            );
        a = max(a, 0.1);
    

	Info << "a-function coefficient: " 
		<< a.weightedAverage(a.mesh().V()).value()
		<<" min: " << min(a).value()
		<<" max: " << max(a).value() << endl;

        volScalarField j3 = 
		(
			beta5() * pow(a, -3) 
			+ pow(a, -1) * beta6() 
			+ a * beta7()
		) * exp(- pow(a,2))
		+ beta8() * pow(a, -4)
		+ beta9() * pow(a, -2)
		+ beta10()
		+ beta11() * pow(a, 2);
        //j3.boundaryField() = 0;
        // deleta force next to wall
        const fvPatchList& patches = j3.mesh().boundary();

        //forAll(patches, patchi)
        //{
            //const fvPatch& curPatch = patches[patchi];
            //if (isType<wallFvPatch>(curPatch))
            //{
                //forAll(curPatch, facei)
                //{
                    //label faceCelli = curPatch.faceCells()[facei];
                    //j3[faceCelli] = 0;
                //}
            //}
        //}
	//-------------------------
	return T_ * j3;
}

tmp<volScalarField> kineticFluidModel::J4() const
{
	dimensionedScalar smallT("smallT", T_.dimensions(), 1e-08);
	//volScalarField a = max( mag(dispersedPhase().U()) * sqrt(3.0 / (8.0 * T_ + smallT)), 0.5);
	volScalarField a = min
            ( 
                mag(dispersedPhase().U()) * sqrt(3.0 / (8.0 * T_ + smallT)),
                20.0
            );
        a = max(a, 0.1);

        volScalarField j4 = 
		(
			beta12() * pow(a, -1) 
			+ a * beta13()
			+ beta14() * pow(a, 3) 
		) * exp(- pow(a,2))
		+ beta15() * pow(a, -2)
		+ beta16()
		+ beta17() * pow(a, 2);
        //j4.boundaryField() = 0;
        // deleta force next to wall
        const fvPatchList& patches = j4.mesh().boundary();

        //forAll(patches, patchi)
        //{
            //const fvPatch& curPatch = patches[patchi];
            //if (isType<wallFvPatch>(curPatch))
            //{
                //forAll(curPatch, facei)
                //{
                    //label faceCelli = curPatch.faceCells()[facei];
                    //j4[faceCelli] = 0;
                //}
            //}
        //}
	//-------------------------

	return T_ * j4;
}

tmp<volVectorField> kineticFluidModel::deltaG()
{
	volScalarField alpha = dispersedPhase();
	volScalarField tau = tau_ -> total();
	volScalarField cd = - fluid_.dragCoeff() 
		/ dispersedPhase().rho() * dispersedPhase();
	volScalarField R = 1.0 / (1.0 + E1_ + 2 * E2_);

	dimensionedScalar epsSmall("epsSmall", epsilon_.dimensions(), 1e-08);
	dimensionedScalar TSmall("TSmall", T_.dimensions(), 1e-08);
	dimensionedScalar cdSmall("cdSmall", cd.dimensions(), 1e-08);
	//dimensionedScalar tauSmall("tauSmall", tau.dimensions(), 1e-03);
	cd = cd + cdSmall;
	//tau = tau + tauSmall;
	

	volScalarField eta1 = 1.0 - 4.0 * tau * cd 
		/ (1.0 - 6.0 * tau * cd + 8.0 * pow(cd, 2) * pow(tau, 2)) + 1e-08;

	volScalarField eta2 = 3 * tau * epsilon_ * (1.0 - 2.0 * tau * cd) / 
		( T_ * (1.0 - 6.0 * tau * cd + 8.0 * pow(cd, 2) * pow(tau, 2))) + 1e-08;

	volScalarField g_cd = eta1 * 2.0 * tau / (1.0 - 2 * tau * cd) * (1.0 - R * E1_ * eta1)
		+ eta2 * 4.0 * tau / (1.0 - 4 * tau * cd) * (1.0 - R * E1_ * eta2)
		+ E2_ * (2.0 + 3.0 / 2.0 * R) *
		( 
			1.0 / cd + 4.0 * tau / (1.0 - 4 * tau * cd) 
		);
	Info << "g_cd coefficient: " 
		<< g_cd.weightedAverage(g_cd.mesh().V()).value()
		<<" min: " << min(g_cd).value()
		<<" max: " << max(g_cd).value() << endl;

	volScalarField g_tau = eta1 * 2.0 * cd / (1.0 - 2 * tau * cd) * (1.0 - R * E1_ * eta1)
		+ eta2 / tau / (1.0 - 4 * tau * cd) * (1.0 - R * E1_ * eta2)
		+ E2_ * (2.0 - 3.0 / 2.0 * R) *
		( 
			1.0 / tau + 4.0 * cd / (1.0 - 4 * tau * cd) 
		);
	Info << "g_tau coefficient: " 
		<< g_tau.weightedAverage(g_tau.mesh().V()).value()
		<<" min: " << min(g_tau).value()
		<<" max: " << max(g_tau).value() << endl;

	volScalarField g_epsilon = eta1 * (1.0 + E1_) / (epsilon_ + epsSmall);
	Info << "g_epsilon coefficient: " 
		<< g_epsilon.weightedAverage(g_epsilon.mesh().V()).value()
		<<" min: " << min(g_epsilon).value()
		<<" max: " << max(g_epsilon).value() << endl;

	volScalarField g_T = - 1.0 / T_ *
		(
			1.0 - E1_ * R - 1.5 * E1_
		)
		+ 3.0 * E2_ / T_;
	Info << "g_T coefficient: " 
		<< g_T.weightedAverage(g_T.mesh().V()).value()
		<<" min: " << min(g_T).value()
		<<" max: " << max(g_T).value() << endl;

	volScalarField g_alpha = (E1_ + 2.0 * E2_) / (alpha + SMALL);
	Info << "g_alpha coefficient: " 
		<< g_alpha.weightedAverage(g_alpha.mesh().V()).value()
		<<" min: " << min(g_alpha).value()
		<<" max: " << max(g_alpha).value() << endl;

        // deleta force next to wall
        const fvPatchList& patches = E2_.mesh().boundary();


        volVectorField g1("g_eps", g_epsilon * fvc::grad(epsilon_));
        volVectorField g2("g_alpha", g_alpha * fvc::grad(alpha));
        volVectorField g3("g_T", g_T * fvc::grad(T_));
        volVectorField g4("g_tau", g_tau * fvc::grad(tau));
        volVectorField g5("g_cd", g_cd * fvc::grad(cd));

	Info << "g-eps: " 
		<< g1.weightedAverage(g_alpha.mesh().V()).value()
		<<" min: " << min(g1).value()
		<<" max: " << max(g1).value() << endl;
	Info << "g-alpha: " 
		<< g2.weightedAverage(g_alpha.mesh().V()).value()
		<<" min: " << min(g2).value()
		<<" max: " << max(g2).value() << endl;
	Info << "g-k: " 
		<< g3.weightedAverage(g_alpha.mesh().V()).value()
		<<" min: " << min(g3).value()
		<<" max: " << max(g3).value() << endl;
	Info << "g-tau: " 
		<< g4.weightedAverage(g_alpha.mesh().V()).value()
		<<" min: " << min(g4).value()
		<<" max: " << max(g4).value() << endl;
	Info << "g-cd: " 
		<< g5.weightedAverage(g_alpha.mesh().V()).value()
		<<" min: " << min(g5).value()
		<<" max: " << max(g5).value() << endl;
        //g1.write();
        //g2.write();
        //g3.write();
        //g4.write();
        //g5.write();

        //g1.boundaryField() = vector(0,0,0);
        //g2.boundaryField() = vector(0,0,0);
        //g3.boundaryField() = vector(0,0,0);
        //g4.boundaryField() = vector(0,0,0);
        //g5.boundaryField() = vector(0,0,0);
        //forAll(patches, patchi)
        //{
                //const fvPatch& curPatch = patches[patchi];
                //if (isType<wallFvPatch>(curPatch))
                //{
                        //forAll(curPatch, facei)
                        //{
                                //label faceCelli = curPatch.faceCells()[facei];
                                //x[faceCelli] = vector(0, 0, 0);
                        //}
                //}
        //}

	return  (g1 + g2 + g3 + g4 + g5);
            //g_alpha * fvc::grad(alpha)
                //+ g_epsilon * fvc::grad(epsilon_)
		//+ g_T * fvc::grad(T_);
                //+ g_tau * fvc::grad(tau);
                //+ g_cd * fvc::grad(cd); 
}


tmp<volVectorField> kineticFluidModel::F1(surfaceScalarField& phi) const
{
	const volVectorField& U = dispersedPhase().U();
	dimensionedScalar smallU("smallU", U.dimensions(), 1e-03);
	const volScalarField rad = g0();
	volScalarField j3 = J3();
	volScalarField j1 = J1();
	volScalarField R = 1.0 / (1.0 + E1_ + 2 * E2_);
	Info << "R: " 
		<< R.weightedAverage(R.mesh().V()).value()
		<<" min: " << min(R).value()
		<<" max: " << max(R).value() << endl;
	Info << "g0: " 
		<< rad.weightedAverage(rad.mesh().V()).value()
		<<" min: " << min(rad).value()
		<<" max: " << max(rad).value() << endl;
	Info << "j3: " 
		<< j3.weightedAverage(j3.mesh().V()).value()
		<<" min: " << min(j3).value()
		<<" max: " << max(j3).value() << endl;
	Info << "j1: " 
		<< j1.weightedAverage(j1.mesh().V()).value()
		<<" min: " << min(j1).value()
		<<" max: " << max(j1).value() << endl;


	volScalarField x = (deltaG_ & U) *
		(
			dispersedPhase().d() * j3 / R / (mag(U) + smallU)
		) + j1;

	return 6.0 * (1.0 + e_) * pow(R, 2) * pow(dispersedPhase(), 2) * g0() * fvc::grad(x);
}

tmp<volVectorField> kineticFluidModel::F2(surfaceScalarField& phi) const
{
	const volVectorField& U = dispersedPhase().U();
	dimensionedScalar smallU("smallU", U.dimensions(), 1e-08);
	volScalarField j3 = J3();
	volScalarField R = 1.0 / (1.0 + E1_ + 2 * E2_);

	volVectorField x = fvc::div(deltaG_ * j3) * U / (mag(U) + smallU);
	x += fvc::div(phi, deltaG_ * j3) / (mag(U) + smallU);
	return 6.0 * (1.0 + e_) * pow(R, 3) * dispersedPhase().d() * pow(dispersedPhase(), 2) * g0() * x;
}

tmp<volVectorField> kineticFluidModel::F3(surfaceScalarField& phi) const
{
	const volVectorField& U = dispersedPhase().U();
	dimensionedScalar smallU("smallU", U.dimensions(), 1e-03);
	volScalarField j1 = J1();
	volScalarField j2 = J2();
	volScalarField j3 = J3();
	volScalarField j4 = J4();
	Info << "j2: " 
		<< j2.weightedAverage(j2.mesh().V()).value()
		<<" min: " << min(j2).value()
		<<" max: " << max(j2).value() << endl;
	Info << "j4: " 
		<< j4.weightedAverage(j4.mesh().V()).value()
		<<" min: " << min(j4).value()
		<<" max: " << max(j4).value() << endl;
	volScalarField R = 1.0 / (1.0 + E1_ + 2 * E2_);

	volScalarField x = (U & deltaG_) * (2.0 * j4 - 5.0 * j3) / (mag(U) + smallU) * R * dispersedPhase().d();
	x += 2.0 * j2 - 3.0 * j1;
	return 6.0 * (1.0 + e_) * pow(R, 2) * pow(dispersedPhase(), 2) * g0() 
		* fvc::div(phi, x) * U / pow(mag(U) + smallU, 2);
}


tmp<volScalarField> kineticFluidModel::g0() const
{
	const volScalarField alpha = dispersedPhase();
	return (2.0 - alpha) / (2.0 * pow(1.0 - alpha, 3));
}

// * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //


}// End namespace Foam
// ************************************************************************* //
