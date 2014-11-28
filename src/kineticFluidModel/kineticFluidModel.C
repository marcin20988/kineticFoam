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
    IOobject::NO_READ,
    IOobject::AUTO_WRITE
   ),
   mesh_,
   dimensionedScalar("E1", dimless, 0.0)
  ),
  E2_
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
   dimensionedScalar("E2", pow(dimLength, -2) / pow(dimTime, -2), 0.0)
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
  )
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

void kineticFluidModel::update()
{
  volScalarField k =  dispersedPhase().turbulence().k();
  volScalarField epsilon =  dispersedPhase().turbulence().epsilon();

  volScalarField tau = tau_ -> field();
  // in the derivation drag coefficient is taken in units of s^-1
  // and is strictly negative
  volScalarField cd = - fluid_.dragCoeff() 
      / dispersedPhase().rho() * dispersedPhase();

  Info << "Collisional relaxation time: " 
    << tau.weightedAverage(tau.mesh().V()).value()
    <<" min: " << min(tau).value()
    <<" max: " << max(tau).value() << endl;
  Info << "Drag coefficient: " << cd.weightedAverage(tau.mesh().V()).value()
    <<" min: " << min(cd).value()
    <<" max: " << max(cd).value() << endl;

  E1_ = 
    (
        2.0 * k - 8 * tau * cd * k 
        + 3.0 * tau * epsilon - 6.0 * cd * pow(tau, 2) * epsilon
    )
    /
    (
        16.0 * pow(cd, 2) * pow(tau, 2) * k - 12.0 * cd * tau * k + 2.0 * k
    );

  E2_ = 3.0 * tau * epsilon / (4.0 * pow(k, 2) * (1.0 - 4.0 * tau * cd));
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
};


tmp<volScalarField> kineticFluidModel::pressureCorrection() const
{
  volScalarField k =  dispersedPhase().turbulence().k();

  return (1.0 + E1_ + 10.0 / 3.0 * k * E2_) / (1.0 + E1_ + 2.0 * k * E2_);
};


tmp<fvVectorMatrix> kineticFluidModel::divDevReff(const volVectorField& U)
{
    volScalarField k =  dispersedPhase().turbulence().k();

    volScalarField correctedViscosity = 8.0 / 9.0 
        * a_ * (1.0 / (1.0 + E1_ + 2.0 * k * E2_)) * pow(k, 2);
    
    volScalarField ratio = 0.5 * correctedViscosity 
        / dispersedPhase().turbulence().nuEff();
    //correctedViscosity -= 2 * dispersedPhase().turbulence().nuEff();
    //correctedViscosity *= dispersedPhase().rho();

    Info << "Ratio of viscosity correction: " << ratio.weightedAverage(k.mesh().V()).value()
        <<" min: " << min(ratio).value()
        <<" max: " << max(ratio).value() << endl;
    //Info << "Corrected viscosity: " << correctedViscosity.weightedAverage(k.mesh().V()).value()
        //<<" min: " << min(correctedViscosity).value()
        //<<" max: " << max(correctedViscosity).value() << endl;

    return
        (
            - fvm::laplacian(correctedViscosity, U)
            - fvc::div(correctedViscosity*dev(T(fvc::grad(U))))
        );
}


tmp<volScalarField> kineticFluidModel::tau() const
{
    return tau_ -> field();
}
// * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //


}// End namespace Foam
// ************************************************************************* //
