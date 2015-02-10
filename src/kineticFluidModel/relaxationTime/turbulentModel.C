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

#include "turbulentModel.H"
#include "addToRunTimeSelectionTable.H"
#include "PhaseIncompressibleTurbulenceModel.H"

namespace Foam
{
namespace relaxationTimes
{
defineTypeNameAndDebug(turbulentModel, 0);
addToRunTimeSelectionTable(relaxationTime, turbulentModel, dictionary);
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

turbulentModel::turbulentModel
(
  const twoPhaseSystem& fluid,
  const dictionary kineticDict,
  const word dispersedPhaseName,
  const phaseModel& dispersedPhase,
  const kineticFluidModel& KM
):
  relaxationTime(fluid, kineticDict, dispersedPhaseName, dispersedPhase, KM),
  minTau_("minTau", dimTime, kineticDict.lookup("minTau")),
  maxTau_("maxTau", dimTime, kineticDict.lookup("maxTau"))
{
};

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //


tmp<volScalarField> turbulentModel::field()
{
  volScalarField nu = dispersedPhase_.turbulence().nuEff();
  volScalarField nut = dispersedPhase_.turbulence().nut() 
      + dimensionedScalar("safeNu", nu.dimensions(), 1e-10);
  volScalarField k =  dispersedPhase_.turbulence().k();
  //volScalarField epsilon =  dispersedPhase_.turbulence().epsilon();

  //volScalarField approximation = max(1.0, nut - nu);
  //Info << "Approximation: " 
    //<< approximation.weightedAverage(k.mesh().V()).value()
    //<<" min: " << min(approximation).value()
    //<<" max: " << max(approximation).value() << endl;

  return max
    (
      27.0 * nu / ( 8.0 * k)
      , minTau_
    );
}

tmp<volScalarField> turbulentModel::total()
{
  volScalarField nu = dispersedPhase_.turbulence().nut();
  volScalarField k = KM_.temp(); //dispersedPhase_.turbulence().k();

  volScalarField A = 
      min
      (
	      27.0 * nu / ( 8.0 * k)
	      * (1.0 + KM_.E1() + 2.0 * KM_.E2()),
	      maxTau_
      );

  tau_ = max(A, minTau_);
  tau_.correctBoundaryConditions();

  return tau_;
}


tmp<volScalarField> turbulentModel::laminar()
{
    return total();
}

tmp<volScalarField> turbulentModel::turbulent()
{
    return total();
}
// * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //


} // end namespace relaxationTimes
} // end namespace Foam
// ************************************************************************* //
