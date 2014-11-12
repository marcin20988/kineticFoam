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

#include "relaxationTime.H"

namespace Foam
{
defineTypeNameAndDebug(relaxationTime, 0);
defineRunTimeSelectionTable(relaxationTime, dictionary);
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<relaxationTime>
relaxationTime::New(const twoPhaseSystem& fluid, const dictionary kineticDict, const word dispersedPhaseName)
{
  word relaxationTimeType(kineticDict.lookup("relaxationTime"));
  Info << "Selecting relaxation time model: " << relaxationTimeType << endl;
  dictionaryConstructorTable::iterator cstrIter = 
    dictionaryConstructorTablePtr_->find(relaxationTimeType);

  if (cstrIter == dictionaryConstructorTablePtr_->end())
  {
    FatalErrorIn("relaxationTime::New")
      << "Unknown relaxationTimeType type "
      << relaxationTimeType << endl << endl
      << "Valid relaxation time types are : " << endl
      << dictionaryConstructorTablePtr_->sortedToc()
      << exit(FatalError);
  }
  return cstrIter()
    (
     fluid,
     kineticDict.subDict(relaxationTimeType+"RelaxationTime"),
     dispersedPhaseName
    );

}


relaxationTime::relaxationTime
(
 const twoPhaseSystem& fluid,
 const dictionary kineticDict,
 const word dispersedPhaseName
)
  //fluid_(fluid),
  //dispersedPhaseName_(dispersedPhaseName)
{
};
// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //
relaxationTime::~relaxationTime(){};


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

void relaxationTime::operator=(const relaxationTime& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorIn("Foam::kineticFluidModel::relaxationTime::operator=(const Foam::kineticFluidModel::relaxationTime&)")
            << "Attempted assignment to self"
            << abort(FatalError);
    }
}

// * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //


} // end namespace Foam
// ************************************************************************* //
