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
#include "PhaseIncompressibleTurbulenceModel.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
  class relaxationTime;

  /*---------------------------------------------------------------------------*\
    Class kineticFluidModel Declaration
    \*---------------------------------------------------------------------------*/

  class kineticFluidModel
  {
    // Private data
    const twoPhaseSystem& fluid_;
    const fvMesh& mesh_;


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

    volScalarField E1_;
    volScalarField E2_;
    volScalarField a_;
    volScalarField T_;

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

    void update(const volScalarField& T);

    tmp<volScalarField> pressureCorrection() const;

    tmp<volScalarField> tau() const;

    tmp<volScalarField> E1() const
    {
        return E1_;
    }

    tmp<volScalarField> E2() const
    {
        return E2_;
    }

    tmp<volScalarField> a() const
    {
        return a_;
    }

    tmp<volScalarField> temp() const
    {
        return T_;
    }

    tmp<fvVectorMatrix> divDevReff(const volVectorField& U);
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
