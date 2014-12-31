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

#include "fullModel.H"
#include "addToRunTimeSelectionTable.H"
#include "PhaseIncompressibleTurbulenceModel.H"

namespace Foam
{
namespace relaxationTimes
{
defineTypeNameAndDebug(fullModel, 0);
addToRunTimeSelectionTable(relaxationTime, fullModel, dictionary);
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

fullModel::fullModel
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


tmp<volScalarField> fullModel::turbulent()
{
  volScalarField nu = dispersedPhase_.turbulence().nut();
  volScalarField k = KM_.temp(); //dispersedPhase_.turbulence().k();
  volScalarField cd = - fluid_.dragCoeff() 
      / KM_.dispersedPhase().rho() * KM_.dispersedPhase() * 0;

  volScalarField A = 
      min
      (
	      27.0 * nu / ( 8.0 * k)
	      * (1.0 + KM_.E1() + 2.0 * k * KM_.E2()),
	      maxTau_
      );

  return max
    (
        A
      , minTau_
    );
}

tmp<volScalarField> fullModel::laminar()
{
  volScalarField nu = dispersedPhase_.nu();
  volScalarField k = KM_.temp(); //dispersedPhase_.turbulence().k();
  volScalarField cd = - fluid_.dragCoeff() 
      / KM_.dispersedPhase().rho() * KM_.dispersedPhase();

  volScalarField A = 
      min
      (
	      27.0 * nu / ( 8.0 * k)
	      * (1.0 + KM_.E1() + 2.0 * k * KM_.E2()),
	      maxTau_
      );

  return max
    (
        A
      , minTau_
    );
}

tmp<volScalarField> fullModel::total()
{
  Info << "Using effective viscosity to calculate relaxation time" << endl;
  volScalarField nu = dispersedPhase_.turbulence().nuEff();
  volScalarField k = KM_.temp(); //dispersedPhase_.turbulence().k();
  volScalarField cd = - fluid_.dragCoeff() 
      / KM_.dispersedPhase().rho() * KM_.dispersedPhase();

  volScalarField A = 
      min
      (
	      27.0 * nu / ( 16.0 * k)
	      * (1.0 + KM_.E1() + 2.0 * k * KM_.E2()),
	      maxTau_
      );

  tau_ = max(A / (1 + A * cd), minTau_);
  tau_.correctBoundaryConditions();

  return tau_;
}
// * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //


} // end namespace relaxationTimes
} // end namespace Foam
// ************************************************************************* //
