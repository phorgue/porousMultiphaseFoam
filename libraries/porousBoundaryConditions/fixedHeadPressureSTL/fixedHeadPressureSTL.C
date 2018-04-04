/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2016 OpenFOAM Foundation
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

#include "fixedHeadPressureSTL.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fixedHeadPressureSTL::
fixedHeadPressureSTL
(
    const fvPatch& h,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(h, iF),
    STLname_("")
{}


Foam::fixedHeadPressureSTL::
fixedHeadPressureSTL
(
    const fvPatch& h,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(h, iF, dict, false),
    STLname_(dict.lookup("file"))
{}


Foam::fixedHeadPressureSTL::
fixedHeadPressureSTL
(
    const fixedHeadPressureSTL& ptf,
    const fvPatch& h,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, h, iF, mapper),
    STLname_(ptf.STLname_)
{}


Foam::fixedHeadPressureSTL::
fixedHeadPressureSTL
(
    const fixedHeadPressureSTL& ptf
)
:
    fixedValueFvPatchScalarField(ptf),
    STLname_(ptf.STLname_)
{}


Foam::fixedHeadPressureSTL::
fixedHeadPressureSTL
(
    const fixedHeadPressureSTL& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF),
    STLname_(ptf.STLname_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fixedHeadPressureSTL::updateCoeffs()
{
    
    if (updated())
    {
        return;
    }
    word STLfile_ = "constant/triSurface/" + STLname_;
    
   triSurfaceMesh potentialSTL(
        IOobject(
            STLfile_,
            this->db()
        )
        );
    pointField pPoints = potentialSTL.points();

    const vectorField& fp = patch().patch().faceCentres();
    scalarField results(patch().patch().faceCentres().size());    
    forAll(fp,facei)
    {
        
        scalar xy_distance = GREAT;
        label id_point = -1;
        forAll(pPoints,pointi)
        {
            scalar tmp_dist = sqrt(pow(pPoints[pointi].x()-fp[facei].x(),2)+pow(pPoints[pointi].y()-fp[facei].y(),2));
            if (tmp_dist < xy_distance)
            {
                xy_distance = tmp_dist;
                id_point = pointi;
            }
        }
        results[facei] = pPoints[id_point].z() - fp[facei].z();
    }
    operator== (results);
    fixedValueFvPatchScalarField::updateCoeffs();
}

void Foam::fixedHeadPressureSTL::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        fixedHeadPressureSTL
    );
}

// ************************************************************************* //
