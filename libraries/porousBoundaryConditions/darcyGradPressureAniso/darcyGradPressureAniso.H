/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
  \\    /   O peration     |
  \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::darcyGradPressureAniso

Description
    Set the pressure gradient boundary condition for a given Darcy total velocity
    Permeability is anisotropic (K volTensorField)

SourceFiles
    darcyGradPressureAniso.C

\*---------------------------------------------------------------------------*/

#ifndef darcyGradPressureAniso_H
#define darcyGradPressureAniso_H

#include "fixedGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
    Class darcyGradPressureFvPatch Declaration
    \*---------------------------------------------------------------------------*/

class darcyGradPressureAniso
:
    public fixedGradientFvPatchScalarField
{
    // Private data

    //- Name of the velocity and permeability fields used to calculate the tensorial wall BC in darcy's law
    word MfName_;
    word MbfName_;
    word phiName_;
    word UName_;
    word LfName_;
    word gradpcName_;		
	
public:

    //- Runtime type information
    TypeName("darcyGradPressureAniso");


    // Constructors

    //- Construct from patch and internal field
    darcyGradPressureAniso
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&
    );

    //- Construct from patch, internal field and dictionary
    darcyGradPressureAniso
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const dictionary&
    );

    //- Construct by mapping given
    //  darcyGradPressureAniso onto a new patch
    darcyGradPressureAniso
    (
        const darcyGradPressureAniso&,
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const fvPatchFieldMapper&
    );

    //- Construct as copy
    darcyGradPressureAniso
    (
        const darcyGradPressureAniso&
    );

    //- Construct and return a clone
    virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
                (
                    new darcyGradPressureAniso(*this)
                );
        }

    //- Construct as copy setting internal field reference
    darcyGradPressureAniso
    (
        const darcyGradPressureAniso&,
        const DimensionedField<scalar, volMesh>&
    );

    //- Construct and return a clone setting internal field reference
    virtual tmp<fvPatchScalarField> clone
    (
        const DimensionedField<scalar, volMesh>& iF
    ) const
        {
            return tmp<fvPatchScalarField>
                (
                    new darcyGradPressureAniso(*this, iF)
                );
        }


    // Member functions

    //- Update the coefficients associated with the patch field
    virtual void updateCoeffs();

    //- Write
    virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
