/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
    Foam::multiscalarMixture
Description
    Multi-component mixture for porous transport.
SourceFiles
    multiscalarMixture.C
\*---------------------------------------------------------------------------*/

#ifndef multiscalarMixture_H
#define multiscalarMixture_H

#include "basicMixture.H"

#include "volFields.H"
#include "PtrList.H"
#include "speciesTable.H"

#include "dispersionModel.H"
#include "sourceEventFile.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class multiscalarMixture Declaration
\*---------------------------------------------------------------------------*/

class multiscalarMixture
:
    public basicMixture
{

public:

    //- Run time type information
    TypeName("multiscalarMixture");


    //- The base class of the mixture
    typedef multiscalarMixture basicMixtureType;


    // Constructors

        //- Construct from dictionary, species names, mesh and phase name
        template<typename EpsType>
        inline multiscalarMixture
        (
            const dictionary& dict,
            const wordList& speciesNames,
            const fvMesh& mesh,
            const word& phaseName,
            const EpsType& eps,
            List<sourceEventFile*>* sourceEventFileRegistry = nullptr,
            const word& sourceEventDtFieldNameOverride = "",
            const dimensionSet& sourceTermDimFactor = dimless
        );

        virtual ~multiscalarMixture()
        {}

    // Member Functions

        //- Return the table of species
        inline const speciesTable& species() const;

        //- Does the mixture include this specie?
        inline bool contains(const word& specieName) const;

        //- Return the mass-fraction fields
        inline PtrList<volScalarField>& Y();

        //- Return the const mass-fraction fields
        inline const PtrList<volScalarField>& Y() const;

        //- Return the mass-fraction field for a specie given by index
        inline volScalarField& Y(const label i);

        //- Return the const mass-fraction field for a specie given by index
        inline const volScalarField& Y(const label i) const;

        //- Return the mass-fraction field for a specie given by name
        inline volScalarField& Y(const word& specieName);

        //- Return the const mass-fraction field for a specie given by name
        const volScalarField& Y(const word& specieName) const;

        //- Retard coefficient [-]
        inline const volScalarField& R(const label speciesi) const;

        //- Decay coefficient [1/s]
        inline const dimensionedScalar& lambda(const label speciesi) const;

        //- Update effective diffusivities and source terms
        template<typename ThetaType>
        inline void correct
        (
            const volVectorField& U,
            const ThetaType& theta
        );

        template<typename SaturationType, typename EpsType>
        inline void correct
        (
            const volVectorField& U,
            const SaturationType& saturation,
            const EpsType& eps
        );


        //- Effective diffusivity [m^2/s]
        inline const volTensorField& Deff(const label speciesi) const;


        //- Source term [Y(speciesi).dimensions()*sourceTermDimFactor/s]
        inline const volScalarField& sourceTerm(const label speciesi) const;

private:
    speciesTable species_;
    PtrList<volScalarField> Y_;
    PtrList<volScalarField> R_;
    List<dimensionedScalar> lambdas_;
    PtrList<dispersionModel> dispersionModels_;
    PtrList<sourceEventFile> sourceEvents_;
    PtrList<volScalarField> sourceTerms_;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "multiscalarMixtureI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
