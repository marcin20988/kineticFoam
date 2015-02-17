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
#include "wallFvPatch.H"

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
    const dimensionedVector g_;


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
    volScalarField R_;
    volScalarField a_;
    volScalarField T_;
    volScalarField epsilon_;
    volVectorField deltaG_;
    volVectorField F_total_;
    volScalarField sp_total_;
    volScalarField NuRatio_;
    // coefficient of restitution
    scalar e_;
    // maximum value of total collisional force component
    scalar maxF_;
    // 'inverse of turbulence intensity' function bounds:
    scalar aMin_;
    scalar aMax_;

    // just to play with; a scalar that multiplies total force
    scalar scaleF_;

    //in deltaG should cd and tau dependent forces be included?
    Switch cdTauForces_;

    // should stokes drag be used instead of drag coefficient from drag model?
    Switch useStokesDrag_;

    // use corrected viscosity? If not turbulent viscosity is used
    Switch useViscosityCorrection_;

    // smooth out collisional force
    Switch forceSmoothing_;

    // delete collisional force on the wall
    Switch wallTreatment_;

    // use deltaG_ forces?
    Switch useG_;

    // development length is a length from inlet over which forces are scaled
    // down to prevent nasty near-inlet effect
    Switch developmentLength_;

    // additional dumping for small forces to prevent oscilations
    Switch dumpSmallF_;

    // additional dumping for small forces to prevent oscilations
    Switch analyticalLimit_;

    //use "turbulence dispersion"?
    Switch dp_;

    // z-coordinate of begining of the development length
    scalar developmentL1_;
    // z-coordinate of end of the development length
    scalar developmentL2_;
    // initial scaling factor on the development length
    scalar developmentScale_;


  public:

    // Static data members

    //- Static data staticData


    // Constructors

    kineticFluidModel(const twoPhaseSystem& fluid, const dimensionedVector& g);

    //- Construct from components

    //- Construct from Istream

    //- Construct as copy


    // Selectors

    //- Select null constructed
    static autoPtr<kineticFluidModel> New(const twoPhaseSystem& fluid, const dimensionedVector g);


    //- Destructor
    ~kineticFluidModel(){};


    // Member Functions
    const phaseModel& dispersedPhase() const;
    const phaseModel& continuousPhase() const;

    void update(const volScalarField& T, const volScalarField& epsilon, const int check = -1);

    tmp<volScalarField> pressureCorrection() const;

    tmp<volScalarField> tauLaminar();

    tmp<volScalarField> tauTurbulent();

    tmp<volScalarField> tauTotal();

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

    tmp<volScalarField> epsilon() const
    {
        return epsilon_;
    }

    tmp<fvVectorMatrix> divDevReff(const volVectorField& U);

    //list of beta functions for J integrals
    tmp<volScalarField> beta(scalar c1, scalar c2, scalar c3) const;
    tmp<volScalarField> beta1() const;
    tmp<volScalarField> beta2() const;
    tmp<volScalarField> beta3() const;
    tmp<volScalarField> beta4() const;
    tmp<volScalarField> beta5() const;
    tmp<volScalarField> beta6() const;
    tmp<volScalarField> beta7() const;
    tmp<volScalarField> beta8() const;
    tmp<volScalarField> beta9() const;
    tmp<volScalarField> beta10() const;
    tmp<volScalarField> beta11() const;
    tmp<volScalarField> beta12() const;
    tmp<volScalarField> beta13() const;
    tmp<volScalarField> beta14() const;
    tmp<volScalarField> beta15() const;
    tmp<volScalarField> beta16() const;
    tmp<volScalarField> beta17() const;

    // J_i functions
    tmp<volScalarField> J1() const;
    tmp<volScalarField> J2() const;
    tmp<volScalarField> J3() const;
    tmp<volScalarField> J4() const;

    // delta g vector
    tmp<volVectorField> deltaG();
	
    // radial distribution function at contact
    tmp<volScalarField> g0() const;

    scalar maxF() const
    {
        return maxF_;
    }
    
    Switch dp() const
    {
        return dp_;
    }

    Switch developmentLength() const
    {
        return developmentLength_;
    }

    tmp<volScalarField> dampFunction()
    {
        dimensionedScalar dL("dL", dimLength, 1.0f);
        volScalarField dampingFunction
            (
                "dampingFunction",
                developmentScale_ 
                + (1.0 - developmentScale_ )
                * pow(mesh_.C().component(2) / dL - developmentL1_, 3)
                / pow(developmentL2_ - developmentL1_, 3)
            );
        return min(dampingFunction, scalar(1.0));
    }

    // this given the force equal to gradient (J1 + u0 dot Delta g R sigma J_3 / u0)
    tmp<volVectorField> F1(surfaceScalarField& phi) const;
    tmp<volVectorField> F2(surfaceScalarField& phi) const;
    tmp<volVectorField> F3(surfaceScalarField& phi) const;
    tmp<volVectorField> F4(surfaceScalarField& phi) const;
    tmp<volVectorField> F5(surfaceScalarField& phi) const;
    tmp<volVectorField> F6(surfaceScalarField& phi) const;
    tmp<volScalarField> F4Sp(surfaceScalarField& phi) const;
    tmp<volScalarField> F5Sp(surfaceScalarField& phi) const;
    tmp<volScalarField> F6Sp(surfaceScalarField& phi) const;
    volVectorField& collisionalF(surfaceScalarField& phi);
    volScalarField& collisionalSp(surfaceScalarField& phi);

    tmp<volVectorField> velocity() const;

    // Access

    // Check

    // Edit

    // Write
    void print();


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
