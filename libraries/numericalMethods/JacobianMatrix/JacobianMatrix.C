/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "JacobianMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

JacobianMatrix::~JacobianMatrix()
{}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


JacobianMatrix::JacobianMatrix
(
    const volScalarField& vf
)
    :
    mesh_(vf.mesh()),
    matrix_(vf, vf.dimensions()),
    lowerListM_(vf.size(),labelList()),
    upperListM_(vf.size(),labelList()),
    lowerListF_(vf.size(),labelList()),
    upperListF_(vf.size(),labelList())
{
    constructAddressing();
}

// * * * * * * * * * * * * * * * Private Members * * * * * * * * * * * * * * //

void JacobianMatrix::constructAddressing()
{
    const lduAddressing& matrixAddressing = matrix_.lduAddr();

    forAll(matrixAddressing.upperAddr(),iter)
    {
        lowerListF_[matrixAddressing.lowerAddr()[iter]].append(matrixAddressing.upperAddr()[iter]);
        lowerListM_[matrixAddressing.lowerAddr()[iter]].append(iter);
    }
    forAll(matrixAddressing.lowerAddr(),iter)
    {
        upperListF_[matrixAddressing.upperAddr()[iter]].append(matrixAddressing.lowerAddr()[iter]);
        upperListM_[matrixAddressing.upperAddr()[iter]].append(iter);
    }
}

void JacobianMatrix::storeColumn(const volScalarField& dF, label celli)
{
    matrix_.diag()[celli] = dF[celli];
    forAll(lowerListM_[celli],iter)
    {
        matrix_.lower()[lowerListM_[celli][iter]] = dF[lowerListF_[celli][iter]];
    }
    forAll(upperListM_[celli],iter)
    {
        matrix_.upper()[upperListM_[celli][iter]] = dF[upperListF_[celli][iter]];
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
