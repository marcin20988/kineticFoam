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
    Foam::kineticFluidModel

Description

SourceFiles
    kineticFluidModelI.H
    kineticFluidModel.C
    kineticFluidModelIO.C

\*---------------------------------------------------------------------------*/

#ifndef kineticFluidModel_H
#define kineticFluidModel_H

#include "twoPhaseSystem.H"
#include "relaxationTime.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

  /*---------------------------------------------------------------------------*\
    Class kineticFluidModel Declaration
    \*---------------------------------------------------------------------------*/

  class kineticFluidModel
  {
    // Private data
    const twoPhaseSystem& fluid_;


    //- Description of data_


    // Private Member Functions

    //- Disallow default bitwise copy construct
    kineticFluidModel(const kineticFluidModel&);

    //- Disallow default bitwise assignment
    void operator=(const kineticFluidModel&);


    protected: 
    const dictionary kineticFluidModelDict_;
    private:
    const word continuousPhaseName_;
    const word dispersedPhaseName_;
    autoPtr<relaxationTime> tau_;

    public:

    // Static data members

    //- Static data staticData


    // Constructors

    kineticFluidModel(const twoPhaseSystem& fluid);

    //- Construct from components

    //- Construct from Istream

    //- Construct as copy


    // Selectors

    //- Select null constructed
    static autoPtr<kineticFluidModel> New(const twoPhaseSystem& fluid);


    //- Destructor
    ~kineticFluidModel(){};


    // Member Functions
    const phaseModel& dispersedPhase() const;

    void update();

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

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
