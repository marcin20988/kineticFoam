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
  const phaseModel& dispersedPhase
):
  relaxationTime(fluid, kineticDict, dispersedPhaseName, dispersedPhase),
  minTau_("minTau", dimTime, kineticDict.lookup("minTau"))
{
};

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //


tmp<volScalarField> turbulentModel::field() const
{
  volScalarField nu = /*dispersedPhase_.nu() + */ dispersedPhase_.turbulence().nuEff();
  volScalarField k =  dispersedPhase_.turbulence().k();
  volScalarField epsilon =  dispersedPhase_.turbulence().epsilon();

  return max
    (
      27.0 * nu * k / ( 8.0 * pow(k,2) - 81.0 * epsilon * nu)
      , minTau_
    );
}

// * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //


} // end namespace relaxationTimes
} // end namespace Foam
// ************************************************************************* //
