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

\*---------------------------------------------------------------------------*/

#include "outputField.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::outputField::outputField
(
    const volScalarField& field,
    const surfaceScalarField& phi,
    const volScalarField& coef1,
    const volScalarField& coef2,
    const volScalarField& coef3,
    bool massBalance,
    bool saturation,
    bool CSVoutput,
    scalar zscale,
    word fileName,
    word dimensions
)
:
    field_(field),
    mesh_(field.mesh()),
    phi_(phi),
    coef1_(coef1),
    coef2_(coef2),
    coef3_(coef3),
    saturation_(saturation),
    massBalance_(massBalance),
    CSVoutput_(CSVoutput),
    zscale_(zscale),
    fileOutput_(fileName)
{
    //- Writing header if CSVoutput is active
    if (CSVoutput) {
        if (saturation) {
            fileOutput_ << "#Time";
        }
        else {
            fileOutput_ << "#Time Total(" << dimensions << ")";
        }
        const fvMesh& mesh = field.mesh();
        forAll(mesh.boundaryMesh(),patchi)
        {
            if (mesh.boundaryMesh()[patchi].type() == "patch")
            {
                fileOutput_ << " flux(" << phi_.boundaryField()[patchi].patch().name() << ")";
            }
        }
        fileOutput_ << endl;
    }
}

template<class Type, template<class> class PatchField, class TypeMesh>
Foam::GeometricField<Type, PatchField, TypeMesh> Foam::outputField::timeInterpolate
(
    const GeometricField<Type, PatchField, TypeMesh>& field,
    const word& timeName,
    scalar ifactor,
    bool writeField
)
{
    GeometricField<Type, PatchField, TypeMesh> ifield
        (
            IOobject
            (
                field.name(),
                timeName,
                field.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            field
        );
    ifield = ifactor*field+(1.0-ifactor)*field.oldTime();
    if (writeField) ifield.write();

    return ifield;
}

void Foam::outputField::write(const word& timeName) {
    if (CSVoutput_) {
        fileOutput_ << timeName << " " << fvc::domainIntegrate(field_ * coef1_ *  coef2_ * coef3_).value()/zscale_;
        if (saturation_) {
            forAll(mesh_.boundaryMesh(), patchi) {
                if (mesh_.boundaryMesh()[patchi].type() == "patch") {
                    fileOutput_ << " " << gSum(phi_.boundaryField()[patchi]);
                }
            }
        }
        else {
            forAll(mesh_.boundaryMesh(), patchi) {
                if (mesh_.boundaryMesh()[patchi].type() == "patch") {
                    fileOutput_ << " " << gSum(field_.boundaryField()[patchi] * phi_.boundaryField()[patchi]);
                }
            }
        }
        fileOutput_ << endl;
    }
}

void Foam::outputField::write(const word& timeName, scalar ifactor) {

    volScalarField fInter = timeInterpolate(field_, timeName, ifactor, true);
    surfaceScalarField phiInter =  timeInterpolate(phi_, timeName, ifactor, true);

    if (CSVoutput_) {
        volScalarField cInter1 = timeInterpolate(coef1_, timeName, ifactor, false);
        volScalarField cInter2 = timeInterpolate(coef2_, timeName, ifactor, false);
        volScalarField cInter3 = timeInterpolate(coef3_, timeName, ifactor, false);
        fileOutput_ << timeName << " " << fvc::domainIntegrate(fInter * cInter1 *  cInter2 * cInter3).value()/zscale_;
        if (saturation_) {
            forAll(mesh_.boundaryMesh(), patchi) {
                if (mesh_.boundaryMesh()[patchi].type() == "patch") {
                    fileOutput_ << " " << gSum(phiInter.boundaryField()[patchi]);
                }
            }
        }
        else {
            forAll(mesh_.boundaryMesh(), patchi) {
                if (mesh_.boundaryMesh()[patchi].type() == "patch") {
                    fileOutput_ << " " << gSum(fInter.boundaryField()[patchi] * phiInter.boundaryField()[patchi]);
                }
            }
        }
        fileOutput_ << endl;
    }
}
