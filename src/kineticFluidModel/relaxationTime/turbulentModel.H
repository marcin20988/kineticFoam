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

Class
    Foam::relaxationTime

Description

SourceFiles
    relaxationTimeI.H
    relaxationTime.C
    relaxationTimeIO.C

\*---------------------------------------------------------------------------*/

#ifndef turbulentModel_H
#define turbulentModel_H

#include "twoPhaseSystem.H"
#include "relaxationTime.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace relaxationTimes
{

/*---------------------------------------------------------------------------*\
                         Class turbulentModel Declaration
\*---------------------------------------------------------------------------*/

class turbulentModel:
  public relaxationTime
{
    // Private data

        //- Description of data_


    // Private Member Functions

        //- Disallow default bitwise copy construct
        //relaxationTime(const relaxationTime&);

        //- Disallow default bitwise assignment
        void operator=(const relaxationTime&);


public:

    // Static data members

        //- Static data staticData


    // Constructors
	TypeName("turbulentModel");

	//- Construct null
	turbulentModel
	  (
	    const twoPhaseSystem& fluid,
	    const dictionary kineticDict,
	    const word dispersedPhaseName,
	    const phaseModel& dispersedPhase
	  );

	//- Construct from components

        //- Construct from Istream

        //- Construct as copy


    // Selectors

        //- Select null constructed


    //- Destructor
    //virtual ~relaxationTime();


    // Member Functions

        // Access

	tmp<volScalarField> field() const;

        // Check

        // Edit

        // Write


    // Member Operators

    // Friend Functions

    // Friend Operators

    // IOstream Operators

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace relaxationTime
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
