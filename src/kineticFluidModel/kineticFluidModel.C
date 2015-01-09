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
  ),
  R_
  (
   IOobject
   (
    "E2",
    mesh_.time().timeName(),
    mesh_,
    IOobject::NO_READ,
    IOobject::NO_WRITE
   ),
   mesh_,
   dimensionedScalar("R",dimless , 0.0)
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
    "deltaG",
    mesh_.time().timeName(),
    mesh_,
    IOobject::MUST_READ,
    IOobject::AUTO_WRITE
   ),
   mesh_
   //dimensionedVector("deltaG", dimless / dimLength, vector(0.0, 0.0, 0.0))
  ),
    F_total_
    (
        IOobject
        (
            "F-total",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
  e_(readScalar(kineticFluidModelDict_.lookup("e"))),
  maxF_(readScalar(kineticFluidModelDict_.lookup("maxF"))),
  aMin_(kineticFluidModelDict_.lookupOrDefault("aMin", 0.1)),
  aMax_(kineticFluidModelDict_.lookupOrDefault("aMax", 20)),
  scaleF_(kineticFluidModelDict_.lookupOrDefault("scaleF", 1.0f)),
  cdTauForces_(kineticFluidModelDict_.lookup("cdTauForces")),
  useStokesDrag_(kineticFluidModelDict_.lookup("useStokesDrag")),
  useViscosityCorrection_(kineticFluidModelDict_.lookup("useViscosityCorrection")),
  forceSmoothing_(kineticFluidModelDict_.lookup("forceSmoothing")),
  wallTreatment_(kineticFluidModelDict_.lookup("wallTreatment")),
  useG_(kineticFluidModelDict_.lookup("useG")),
  developmentLength_(kineticFluidModelDict_.lookup("developmentLength")),
  developmentL1_(readScalar(kineticFluidModelDict_.lookup("developmentLstart"))),
  developmentL2_(readScalar(kineticFluidModelDict_.lookup("developmentLend"))),
  developmentScale_(readScalar(kineticFluidModelDict_.lookup("developmentScale")))
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
  //Info << "Updating kinetic model" << endl;
  T_ = T;
  epsilon_ = epsilon;
  volScalarField k =  T;

  // value of check dictates which tau is to be used
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

  // Stokes drag formula
  volScalarField cd_stokes = - 3.0 * 3.14 
      * fluid_.otherPhase(dispersedPhase()).nu()
      * fluid_.otherPhase(dispersedPhase()).d() 
      / (3.13 / 6.0 * pow(fluid_.otherPhase(dispersedPhase()).d(), 3));



  //Info << "Stokes drag: " << cd_stokes.weightedAverage(tau.mesh().V()).value() << endl;
  //Info << "Drag coefficient: " << cd.weightedAverage(tau.mesh().V()).value()
    //<<" min: " << min(cd).value()
    //<<" max: " << max(cd).value() << endl;
  if(useStokesDrag_)
  {
    //Info << "Using stokes drag" << endl;
    cd = cd_stokes;
  }

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

  R_ = 1.0 + E1_ + E2_;


  a_ = 2.0 * tau / (k *(1.0 - 6.0 * tau * cd));


  deltaG_ = deltaG();
  deltaG_.correctBoundaryConditions();


};


tmp<volScalarField> kineticFluidModel::pressureCorrection() const
{
  volScalarField k =  temp();

  return (1.0 + E1_ + 10.0 / 3.0 * E2_) / (1.0 + E1_ + 2.0 * E2_);
};


tmp<fvVectorMatrix> kineticFluidModel::divDevReff(const volVectorField& U)
{
    volScalarField k =  temp();

    volScalarField correctedViscosity = 0.5 * 8.0 / 9.0 
        * a_ * (1.0 / (1.0 + E1_ + 2.0 * E2_)) * pow(k, 2);

    correctedViscosity += dispersedPhase().nu();

    if(!useViscosityCorrection_)
    {
        correctedViscosity = dispersedPhase().turbulence().nut();
    }

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
                aMax_
            );
        a = max(a, aMin_);

	Info << "a-function coefficient: " 
		<< a.weightedAverage(a.mesh().V()).value()
		<<" min: " << min(a).value()
		<<" max: " << max(a).value() << endl;
        volScalarField j1 = beta1() + pow(a, 2) * beta2();
        return min
            (
                // this is an analytical J1
                -T_ * j1,
                // and this is limit og J1 as k -> 0
                beta(-0.08, -0.68, -0.96) * pow
                (
                    mag(dispersedPhase().U()), 
                    2
                )
            );
}

tmp<volScalarField> kineticFluidModel::J2() const
{
	dimensionedScalar smallT("smallT", T_.dimensions(), 1e-08);
	volScalarField a = min
            ( 
                mag(dispersedPhase().U()) * sqrt(3.0 / (8.0 * T_ + smallT)),
                aMax_
            );
        a = max(a, aMin_);

        volScalarField j2 = beta3() + pow(a, 2) * beta4();

	return - T_ * j2;
}

tmp<volScalarField> kineticFluidModel::J3() const
{
	dimensionedScalar smallT("smallT", T_.dimensions(), 1e-08);
	volScalarField a = min
            ( 
                mag(dispersedPhase().U()) * sqrt(3.0 / (8.0 * T_ + smallT)),
                aMax_
            );
        a = max(a, aMin_);
    
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

	return T_ * j3;
}

tmp<volScalarField> kineticFluidModel::J4() const
{
	dimensionedScalar smallT("smallT", T_.dimensions(), 1e-08);
	volScalarField a = min
            ( 
                mag(dispersedPhase().U()) * sqrt(3.0 / (8.0 * T_ + smallT)),
                aMax_
            );
        a = max(a, aMin_);

        volScalarField j4 = 
		(
			beta12() * pow(a, -1) 
			+ a * beta13()
			+ beta14() * pow(a, 3) 
		) * exp(- pow(a,2))
		+ beta15() * pow(a, -2)
		+ beta16()
		+ beta17() * pow(a, 2);

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
	cd = cd + cdSmall;
	

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
        // if switch is off set this coefficient to zero
        if(!cdTauForces_)
        {
            g_cd = g_cd * 0;
        }

	volScalarField g_tau = eta1 * 2.0 * cd / (1.0 - 2 * tau * cd) 
                * (1.0 - R * E1_ * eta1)
		+ eta2 / tau / (1.0 - 4 * tau * cd) * (1.0 - R * E1_ * eta2)
		+ E2_ * (2.0 - 3.0 / 2.0 * R) *
		( 
			1.0 / tau + 4.0 * cd / (1.0 - 4 * tau * cd) 
		);
        // if switch is off set this coefficient to zero
        if(!cdTauForces_)
        {
            g_tau = 0.0 * g_tau;
        }

	volScalarField g_epsilon = eta1 * (1.0 + E1_);

        volScalarField logEps = 
            log(
                max(epsilon_, epsSmall)
                / dimensionedScalar("eps", epsilon_.dimensions(), 1.0)
            );

	volScalarField g_T = - 
		(
			1.0 - E1_ * R - 1.5 * E1_- 3.0 * E2_
		);
        volScalarField logT = 
            log
            (
                max(T_, TSmall)
                / dimensionedScalar("T", T_.dimensions(), 1.0)
            );

	volScalarField g_alpha = (E1_ + 2.0 * E2_);
        volScalarField logAlpha = log(max(alpha, 1e-08));
        //----------------------------------------------------------------------
	//Info << "g_cd coefficient: " 
		//<< g_cd.weightedAverage(g_cd.mesh().V()).value()
		//<<" min: " << min(g_cd).value()
		//<<" max: " << max(g_cd).value() << endl;
	//Info << "g_tau coefficient: " 
		//<< g_tau.weightedAverage(g_tau.mesh().V()).value()
		//<<" min: " << min(g_tau).value()
		//<<" max: " << max(g_tau).value() << endl;
	//Info << "g_epsilon coefficient: " 
		//<< g_epsilon.weightedAverage(g_epsilon.mesh().V()).value()
		//<<" min: " << min(g_epsilon).value()
		//<<" max: " << max(g_epsilon).value() << endl;
	//Info << "g_T coefficient: " 
		//<< g_T.weightedAverage(g_T.mesh().V()).value()
		//<<" min: " << min(g_T).value()
		//<<" max: " << max(g_T).value() << endl;
	//Info << "g_alpha coefficient: " 
		//<< g_alpha.weightedAverage(g_alpha.mesh().V()).value()
		//<<" min: " << min(g_alpha).value()
		//<<" max: " << max(g_alpha).value() << endl;
        //----------------------------------------------------------------------
        volVectorField g1("g_eps", g_epsilon * fvc::grad(logEps));
        volVectorField g2("g_alpha", g_alpha * fvc::grad(logAlpha));
        volVectorField g3("g_T", g_T * fvc::grad(logT));
        volVectorField g4("g_tau", g_tau * fvc::grad(tau));
        volVectorField g5("g_cd", g_cd * fvc::grad(cd));

	//Info << "g-eps: " 
		//<< g1.weightedAverage(g_alpha.mesh().V()).value()
		//<<" min: " << min(g1).value()
		//<<" max: " << max(g1).value() << endl;
	//Info << "g-alpha: " 
		//<< g2.weightedAverage(g_alpha.mesh().V()).value()
		//<<" min: " << min(g2).value()
		//<<" max: " << max(g2).value() << endl;
	//Info << "g-k: " 
		//<< g3.weightedAverage(g_alpha.mesh().V()).value()
		//<<" min: " << min(g3).value()
		//<<" max: " << max(g3).value() << endl;
	//Info << "g-tau: " 
		//<< g4.weightedAverage(g_alpha.mesh().V()).value()
		//<<" min: " << min(g4).value()
		//<<" max: " << max(g4).value() << endl;
	//Info << "g-cd: " 
		//<< g5.weightedAverage(g_alpha.mesh().V()).value()
		//<<" min: " << min(g5).value()
		//<<" max: " << max(g5).value() << endl;

	return  (g1 + g2 + g3 + g4 + g5);
}


tmp<volVectorField> kineticFluidModel::F1(surfaceScalarField& phi) const
{
	const volVectorField& U = dispersedPhase().U();
	dimensionedScalar smallU("smallU", U.dimensions(), 1e-03);

	const volScalarField rad = g0();

	volScalarField j1("J1", J1());

	return fvc::grad
            ( 
                6.0 * (1.0 + e_) * pow(R_, 2) * pow(dispersedPhase(), 2) 
                * g0() * j1
            );
}

tmp<volVectorField> kineticFluidModel::F2(surfaceScalarField& phi) const
{
	const volVectorField& U = dispersedPhase().U();
	dimensionedScalar smallU("smallU", U.dimensions(), 1e-03);

	const volScalarField rad = g0();

	volScalarField j3 = J3();

	volScalarField x
            (
                "x",
                (deltaG_ & U) *
                (
                    dispersedPhase().d() * j3 / R_ / (mag(U) + smallU)
                )
            );

	return fvc::grad
            (
                6.0 * (1.0 + e_) * pow(R_, 2) * pow(dispersedPhase(), 2) 
                * g0() * x
            );
}

tmp<volVectorField> kineticFluidModel::F3(surfaceScalarField& phi) const
{
	const volVectorField& U = dispersedPhase().U();
	dimensionedScalar smallU("smallU", U.dimensions(), 1e-08);
	volScalarField j3("J3", J3());

	//volVectorField x("x", fvc::div(phi, deltaG_ * j3) / (mag(U) + smallU));
	//return 6.0 * (1.0 + e_) * pow(R_, 3) * dispersedPhase().d() 
            //* pow(dispersedPhase(), 2) * g0() * x;
	volVectorField x
            (
                "x", 
                deltaG_ * j3 
                * 6.0 * (1.0 + e_) * pow(R_, 3) * dispersedPhase().d() 
                * pow(dispersedPhase(), 2) * g0()
            );
	return fvc::div(phi, x) / (mag(U) + smallU);
}

tmp<volVectorField> kineticFluidModel::F4(surfaceScalarField& phi) const
{
	const volVectorField& U = dispersedPhase().U();
	dimensionedScalar smallU("smallU", U.dimensions(), 1e-08);
	volScalarField j3("J3", J3());

	//volVectorField x = fvc::div(deltaG_ * j3) * U / (mag(U) + smallU);
	//return 6.0 * (1.0 + e_) * pow(R_, 3) * dispersedPhase().d() 
            //* pow(dispersedPhase(), 2) * g0() * x;
        volVectorField x
            (
                "x",
                deltaG_ * j3
                * 6.0 * (1.0 + e_) * pow(R_, 3) * dispersedPhase().d() 
                * pow(dispersedPhase(), 2) * g0()
            );
        return fvc::div(x) * U / (mag(U) + smallU);
}

tmp<volScalarField> kineticFluidModel::F4Sp(surfaceScalarField& phi) const
{
	const volVectorField& U = dispersedPhase().U();
	dimensionedScalar smallU("smallU", U.dimensions(), 1e-08);
	volScalarField j3("J3", J3());

	//volVectorField x = fvc::div(deltaG_ * j3) * U / (mag(U) + smallU);
	//return 6.0 * (1.0 + e_) * pow(R_, 3) * dispersedPhase().d() 
            //* pow(dispersedPhase(), 2) * g0() * x;
        volVectorField x
            (
                "x",
                deltaG_ * j3
                * 6.0 * (1.0 + e_) * pow(R_, 3) * dispersedPhase().d() 
                * pow(dispersedPhase(), 2) * g0()
            );
        return fvc::div(x) / (mag(U) + smallU);
}


tmp<volVectorField> kineticFluidModel::F5(surfaceScalarField& phi) const
{
	const volVectorField& U = dispersedPhase().U();
	dimensionedScalar smallU("smallU", U.dimensions(), 1e-03);
	volScalarField j3 = J3();
	volScalarField j4 = J4();

	//volScalarField x
            //(
                //"x",
                //(U & deltaG_) * (2.0 * j4 - 5.0 * j3) 
                /// (mag(U) + smallU) * R_ * dispersedPhase().d()
            //);
	//return 6.0 * (1.0 + e_) * pow(R_, 2) * pow(dispersedPhase(), 2) * g0() 
		//* fvc::div(phi, x) * U / pow(mag(U) + smallU, 2);
	volScalarField x
            (
                "x",
                6.0 * (1.0 + e_) * pow(R_, 2) * pow(dispersedPhase(), 2) * g0() 
                * (U & deltaG_) * (2.0 * j4 - 5.0 * j3) 
                / (mag(U) + smallU) * R_ * dispersedPhase().d()
            );
	return (fvc::grad(x) & U) * U / pow(mag(U) + smallU, 2);
}

tmp<volScalarField> kineticFluidModel::F5Sp(surfaceScalarField& phi) const
{
	const volVectorField& U = dispersedPhase().U();
	dimensionedScalar smallU("smallU", U.dimensions(), 1e-03);
	volScalarField j3 = J3();
	volScalarField j4 = J4();

	//volScalarField x
            //(
                //"x",
                //(U & deltaG_) * (2.0 * j4 - 5.0 * j3) 
                /// (mag(U) + smallU) * R_ * dispersedPhase().d()
            //);
	//return 6.0 * (1.0 + e_) * pow(R_, 2) * pow(dispersedPhase(), 2) * g0() 
		//* fvc::div(phi, x) * U / pow(mag(U) + smallU, 2);
	volScalarField x
            (
                "x",
                6.0 * (1.0 + e_) * pow(R_, 2) * pow(dispersedPhase(), 2) * g0() 
                * (U & deltaG_) * (2.0 * j4 - 5.0 * j3) 
                / (mag(U) + smallU) * R_ * dispersedPhase().d()
            );
	return (fvc::grad(x) & U) / pow(mag(U) + smallU, 2);
}

tmp<volVectorField> kineticFluidModel::F6(surfaceScalarField& phi) const
{
	const volVectorField& U = dispersedPhase().U();
	dimensionedScalar smallU("smallU", U.dimensions(), 1e-03);
	volScalarField j1 = J1();
	volScalarField j2 = J2();

	volScalarField x("x", 2.0 * j2 - 3.0 * j1);
	return fvc::grad
            (
                //phi,
                6.0 * (1.0 + e_) * pow(R_, 2) * pow(dispersedPhase(), 2) 
                * g0() *  x
            ) & U * U / pow(mag(U) + smallU, 2);
}

tmp<volScalarField> kineticFluidModel::F6Sp(surfaceScalarField& phi) const
{
	const volVectorField& U = dispersedPhase().U();
	dimensionedScalar smallU("smallU", U.dimensions(), 1e-03);
	volScalarField j1 = J1();
	volScalarField j2 = J2();

	volScalarField x("x", 2.0 * j2 - 3.0 * j1);
	return fvc::grad
            (
                //phi,
                6.0 * (1.0 + e_) * pow(R_, 2) * pow(dispersedPhase(), 2) 
                * g0() *  x
            ) & U / pow(mag(U) + smallU, 2);
}



tmp<volScalarField> kineticFluidModel::g0() const
{
	const volScalarField alpha = dispersedPhase();
	return (2.0 - alpha) / (2.0 * pow(1.0 - alpha, 3));
}

volVectorField& kineticFluidModel::collisionalF(surfaceScalarField& phi)
{
    volVectorField f1("F1", F1(phi));
    volVectorField f2("F2", F2(phi));
    volVectorField f3("F3", F3(phi));
    //volVectorField f4("F4", F4(phi));
    //volVectorField f5("F5", F5(phi));
    //volVectorField f6("F6", F6(phi));

    volScalarField alpha = dispersedPhase();

    // I think F1, F2, F4 needs to be multiplied bu alpha to get phase-sensitive
    // formulation of OF momentum equation
    // other possibility is F3, F5, F6 should be divided by alpha
    // but I think first one is correct

    F_total_ = alpha * f1;
    if(useG_)
    {
        F_total_ += alpha * (f2) + f3;
    }
    F_total_ *= scaleF_;

    F_total_.correctBoundaryConditions();

    forAll(mesh_.C(), celli)
    {
        F_total_[celli].x() = min(F_total_[celli].x(), maxF_);
        F_total_[celli].y() = min(F_total_[celli].y(), maxF_);
        F_total_[celli].z() = min(F_total_[celli].z(), maxF_);

        F_total_[celli].x() = max(F_total_[celli].x(), -maxF_);
        F_total_[celli].y() = max(F_total_[celli].y(), -maxF_);
        F_total_[celli].z() = max(F_total_[celli].z(), -maxF_);
    }

    volVectorField F0("F0", F_total_);

    if(developmentLength_)
    {

        forAll(mesh_.C(), celli)
        {
            scalar coord = mesh_.C()[celli].z();

            if(coord < developmentL2_ && coord > developmentL1_)
            {
                F0[celli] *= 
                    (
                        developmentScale_ 
                        + (1.0 - developmentScale_ )
                        * pow(coord - developmentL1_, 3)
                        / pow(developmentL2_ - developmentL1_, 3)
                    );
                F_total_[celli] *= 
                    (
                        developmentScale_ 
                        + (1.0 - developmentScale_ )
                        * pow(coord - developmentL1_, 3)
                        / pow(developmentL2_ - developmentL1_, 3)
                    );
            }
        }
    }
    
    if(wallTreatment_)
    {
        const fvPatchList& patches = mesh_.boundary();

        forAll(patches, patchi)
        {
                const fvPatch& curPatch = patches[patchi];
                if (isType<wallFvPatch>(curPatch))
                {
                        forAll(curPatch, facei)
                        {
                                label faceCelli = curPatch.faceCells()[facei];
                                F0[faceCelli] = vector(0, 0, 0);
                                F_total_[faceCelli] = vector(0, 0, 0);
                        }
                }
        }
    }

    if(forceSmoothing_)
    {
        forAll(mesh_.C(), celli)
        {
                scalar count = 1.0;
                forAll(mesh_.cellCells()[celli], cellj)
                {
                    label n_id = mesh_.cellCells()[celli][cellj];
                    count += 1.0;
                    F_total_[celli] += F0[n_id];
                }
                F_total_[celli] /= count;
        }
    }

    if(mesh_.time().outputTime())
    {
        f1.write();
        f2.write();
        f3.write();
        //f4.write();
        //f5.write();
        //f6.write();
    }

    F_total_.correctBoundaryConditions();
    return F_total_;
}


tmp<volScalarField> kineticFluidModel::collisionalSp(surfaceScalarField& phi)
{
    volScalarField f6("F6", F6Sp(phi));
    volScalarField f4("F4", F4Sp(phi));
    volScalarField f5("F5", F5Sp(phi));

    volScalarField alpha = dispersedPhase();

    // I think F1, F2, F4 needs to be multiplied bu alpha to get phase-sensitive
    // formulation of OF momentum equation
    // other possibility is F3, F5, F6 should be divided by alpha
    // but I think first one is correct

    volScalarField sp_total_ = alpha * f6;
    if(useG_)
    {
        sp_total_ += alpha * (/*f4 +*/ f5);
    }
    sp_total_ *= scaleF_;

    sp_total_.boundaryField() = 0;

    forAll(mesh_.C(), celli)
    {
        sp_total_[celli] = min(sp_total_[celli], maxF_);
    }

    volScalarField F0("F0", sp_total_);

    if(developmentLength_)
    {

        forAll(mesh_.C(), celli)
        {
            scalar coord = mesh_.C()[celli].z();

            if(coord < developmentL2_ && coord > developmentL1_)
            {
                F0[celli] *= 
                    (
                        developmentScale_ 
                        + (1.0 - developmentScale_ )
                        * pow(coord - developmentL1_, 3)
                        / pow(developmentL2_ - developmentL1_, 3)
                    );
                sp_total_[celli] *= 
                    (
                        developmentScale_ 
                        + (1.0 - developmentScale_ )
                        * pow(coord - developmentL1_, 3)
                        / pow(developmentL2_ - developmentL1_, 3)
                    );
            }
        }
    }
    
    if(wallTreatment_)
    {
        const fvPatchList& patches = mesh_.boundary();

        forAll(patches, patchi)
        {
                const fvPatch& curPatch = patches[patchi];
                if (isType<wallFvPatch>(curPatch))
                {
                        forAll(curPatch, facei)
                        {
                                label faceCelli = curPatch.faceCells()[facei];
                                F0[faceCelli] = 0;
                                sp_total_[faceCelli] = 0;
                        }
                }
        }
    }

    if(forceSmoothing_)
    {
        forAll(mesh_.C(), celli)
        {
                scalar count = 1.0;
                forAll(mesh_.cellCells()[celli], cellj)
                {
                    label n_id = mesh_.cellCells()[celli][cellj];
                    count += 1.0;
                    sp_total_[celli] += F0[n_id];
                }
                sp_total_[celli] /= count;
        }
    }

    if(mesh_.time().outputTime())
    {
        f6.write();
        f4.write();
        f5.write();
    }

    sp_total_.boundaryField() = 0;
    return sp_total_ * 1.0;
}



void kineticFluidModel::print()
{
  volScalarField tau = tau_ -> total();

  Info << "Collisional relaxation time: " 
    << tau.weightedAverage(tau.mesh().V()).value()
    <<" min: " << min(tau).value()
    <<" max: " << max(tau).value() << endl;

  Info << "E1: " << E1_.weightedAverage(tau.mesh().V()).value()
    <<" min: " << min(E1_).value()
    <<" max: " << max(E1_).value() << endl;
  Info << "E2: " << E2_.weightedAverage(tau.mesh().V()).value()
    <<" min: " << min(E2_).value()
    <<" max: " << max(E2_).value() << endl;
  Info << "R: " << R_.weightedAverage(tau.mesh().V()).value()
    <<" min: " << min(R_).value()
    <<" max: " << max(R_).value() << endl;
  
  Info << "a: " << a_.weightedAverage(tau.mesh().V()).value()
    <<" min: " << min(a_).value()
    <<" max: " << max(a_).value() << endl;

  Info << "deltaG: " << deltaG_.weightedAverage(tau.mesh().V()).value()
      <<" min: " << min(deltaG_).value()
      <<" max: " << max(deltaG_).value() << endl;

  Info << "collisional force: " << F_total_.weightedAverage(tau.mesh().V()).value()
      <<" min: " << min(F_total_).value()
      <<" max: " << max(F_total_).value() << endl;
}
// * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //


}// End namespace Foam
// ************************************************************************* //
