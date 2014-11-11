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

#ifndef relaxationTime_H
#define relaxationTime_H

#include "twoPhaseSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticModel
{

/*---------------------------------------------------------------------------*\
                         Class relaxationTime Declaration
\*---------------------------------------------------------------------------*/

class relaxationTime
{
    // Private data

        //- Description of data_


    // Private Member Functions

        //- Disallow default bitwise copy construct
        relaxationTime(const relaxationTime&);

        //- Disallow default bitwise assignment
        void operator=(const relaxationTime&);


public:

    // Static data members

        //- Static data staticData


    // Constructors

        //- Construct null
        relaxationTime();

        //- Construct from components

        //- Construct from Istream

        //- Construct as copy


    // Selectors

        //- Select null constructed
        static autoPtr<relaxationTime> New();


    //- Destructor
    ~relaxationTime();


    // Member Functions

        // Access

        // Check

        // Edit

        // Write


    // Member Operators

    // Friend Functions

    // Friend Operators

    // IOstream Operators

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace kineticModel
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
