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
     dispersedPhase()
   )
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
  volScalarField tau = tau_ -> field();

  Info << "Collisional relaxation time: " << tau.weightedAverage(tau.mesh().V()).value()
    <<" min: " << min(tau).value()
    <<" max: " << max(tau).value() << endl;
};
// * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //


}// End namespace Foam
// ************************************************************************* //
