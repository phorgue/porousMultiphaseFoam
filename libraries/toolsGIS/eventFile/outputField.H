/*---------------------------------------------------------------------------*\
  		  _______  ____    ____  ________
 		 |_   __ \|_   \  /   _||_   __  |
   		   | |__) | |   \/   |    | |_ \_|
   		   |  ___/  | |\  /| |    |  _|
    		  _| |_    _| |_\/_| |_  _| |_
   		 |_____|  |_____||_____||_____|
   	     Copyright (C) Toulouse INP, Pierre Horgue

License
    This file is part of porousMultiphaseFoam, an extension of OpenFOAM
    developed by Pierre Horgue (phorgue@imft.fr) and dedicated to multiphase
    flows through porous media.

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
    Foam::outputField

Description
    Class for outputEvent for time interpolation + mass balance

 SourceFiles
    outputField.C

\*---------------------------------------------------------------------------*/

#ifndef outputField_H
#define outputField_H

#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "OFstream.H"

namespace Foam {

/*---------------------------------------------------------------------------* \
    Class outputEventFile Declaration
\*---------------------------------------------------------------------------*/

class outputField {

private:
    const volScalarField& field_;
    const fvMesh& mesh_;
    const surfaceScalarField& phi_;
    const volScalarField& coef1_;
    const volScalarField& coef2_;
    const volScalarField& coef3_;
    const bool saturation_;
    const bool CSVoutput_;
    const scalar zscale_;
    const word filename_;
    autoPtr<OFstream> fileStream_;
    bool hasSourceTerm_;
    List<word> sourceNames_;
    List<const scalar*> sourceValues_;
    List<scalar> sourceOldValues_;

    template<class Type, template<class> class PatchField, class TypeMesh> GeometricField<Type, PatchField, TypeMesh>
            timeInterpolate(
                    const GeometricField<Type, PatchField, TypeMesh>&,
                    const word& timeName,
                    scalar ifactor,
                    bool writeField = true
                    );

public:
    //- constructor
    outputField(
        const volScalarField& field,
        const surfaceScalarField& phi,
        const volScalarField& coef1,
        const volScalarField& coef2,
        const volScalarField& coef3,
        bool saturation,
        bool CSVoutput,
        const scalar zscale,
        const word& fileName
    );

    //- initialization functions
    void writeHeader();
    void addSourceTerm(const word& name, const scalar& value);

    //- access functions
    const fvMesh &mesh() { return field_.mesh(); };
    const volScalarField &field() { return field_; };

    //- output function
    void write(const scalar& timeValue);
    void write(const word& timeName, const scalar& timeValue, scalar ifactor);

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
