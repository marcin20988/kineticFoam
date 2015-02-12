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

#include "constantModel.H"
#include "addToRunTimeSelectionTable.H"
#include "PhaseIncompressibleTurbulenceModel.H"

namespace Foam
{
namespace relaxationTimes
{
defineTypeNameAndDebug(constantModel, 0);
addToRunTimeSelectionTable(relaxationTime, constantModel, dictionary);
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

constantModel::constantModel
(
  const twoPhaseSystem& fluid,
  const dictionary kineticDict,
  const word dispersedPhaseName,
  const phaseModel& dispersedPhase,
  const kineticFluidModel& KM
):
  relaxationTime(fluid, kineticDict, dispersedPhaseName, dispersedPhase, KM),
  value_("value", dimTime, kineticDict.lookup("tauValue"))
{
};

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //


tmp<volScalarField> constantModel::turbulent()
{
    return total();
}

tmp<volScalarField> constantModel::laminar()
{
    return total();
}

tmp<volScalarField> constantModel::total()
{
 
  volScalarField value
      (
          IOobject
          (
              "dpdt",
              KM_.dispersedPhase().mesh().time().constant(),
              KM_.dispersedPhase().mesh()
          ),
          KM_.dispersedPhase().mesh(),
          value_
          //dimensionedScalar("value", dimTime, value_)
      );

  return value * 1.0;
}
// * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //


} // end namespace relaxationTimes
} // end namespace Foam
// ************************************************************************* //
