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

Namespace
    Foam::dispersionModel

Description
    A namespace for various dispersion implementations.

Class
    Foam::dispersionModel

Description
    An abstract base class for dispersion models.

SourceFiles
    dispersionModel.C
    dispersionModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef dispersionModel_H
#define dispersionModel_H

#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class dispersionModel Declaration
\*---------------------------------------------------------------------------*/

class dispersionModel
{
protected:

    // Protected data
    word name_;
    dictionary transportProperties_;

    //-Reference to the mesh
    const fvMesh& mesh_;

    //- Molecular and effective diffusion
    dimensionedScalar Dm_;
    volTensorField Deff_;

    //- Disallow copy construct
    dispersionModel(const dispersionModel&);

    //- Disallow default bitwise assignment
    void operator=(const dispersionModel&);


public:

    //- Runtime type information
    TypeName("dispersionModel");


    // Declare run-time constructor selection table

    declareRunTimeSelectionTable
    (
        autoPtr,
        dispersionModel,
        dictionary,
        (
            const word& name,
            const dictionary& transportProperties,
            const fvMesh& mesh
        ),
        (name, transportProperties, mesh)
    );


    // Selectors

    //- Return a reference to the selected dispersion model
    static autoPtr<dispersionModel> New
    (
        const word& name,
        const dictionary& transportProperties,
        const fvMesh& mesh
    );


    // Constructors

    //- Construct from components
    dispersionModel
    (
        const word& name,
        const dictionary& transportProperties,
        const fvMesh& mesh
    );

    //- Destructor
    virtual ~dispersionModel()
        {}


    // Member Functions

    //- Return the phase transport properties dictionary
    const dictionary& transportProperties() const
        {
            return transportProperties_;
        }

    //- Return molecular diffusion
    const dimensionedScalar& Dm() const
        {
            return Dm_;
        }

    //- Return effective dispersion
    tmp<volTensorField> Deff() const
        {
            return Deff_;
        }

    //- Correct the dispersion giving water content or porosity field
    virtual void correct(
        const volScalarField& C,
        const volVectorField& U,
        const volScalarField& theta
    ) = 0;

    //- Correct the dispersion giving scalar porosity or water content
    virtual void correct(
        const volScalarField& C,
        const volVectorField& U,
        const dimensionedScalar& eps
    )
        {
            volScalarField theta
            (
                IOobject
                (
                    "theta",
                    C.time().constant(),
                    C.mesh(),
                    IOobject::NO_READ
                ),
                C.mesh(),
                eps
            );

            return correct(C, U, theta);
        }

    //- Correct the dispersion giving both saturation and porosity (as scalars or fields)
    template<typename SaturationType, typename EpsType>
    void correct(
        const volScalarField& C,
        const volVectorField& U,
        const SaturationType& saturation,
        const EpsType& eps
    )
        {
            correct(C, U, eps*saturation);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
