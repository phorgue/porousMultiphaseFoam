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
    Foam::relativePermeabilityModel

Description
    A namespace for various relative permeability implementations.

Class
    Foam::relativePermeabilityModel

Description
    An abstract base class for relative permeability models.

SourceFiles
    relativePermeabilityModel.C
    relativePermeabilityModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef relativePermeabilityModel_H
#define relativePermeabilityModel_H

#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class relativePermeabilityModel Declaration
\*---------------------------------------------------------------------------*/

class relativePermeabilityModel
{
protected:

    // Protected data
    //- name of the saturation field
    word Sname_;

    // Reference to the dictionary
    const dictionary& transportProperties_;

    //- relative permeability and derivative for each phase
    volScalarField kra_,krb_;
    volScalarField dkradS_,dkrbdS_;

    //- Effective saturation of phase b
    volScalarField Se_;

    //- Disallow copy construct
    relativePermeabilityModel(const relativePermeabilityModel&);

    //- Disallow default bitwise assignment
    void operator=(const relativePermeabilityModel&);


public:

    //- Runtime type information
    TypeName("relativePermeabilityModel");


    // Declare run-time constructor selection table

    declareRunTimeSelectionTable
    (
        autoPtr,
        relativePermeabilityModel,
        dictionary,
        (
            const fvMesh& mesh,
            const dictionary& transportProperties,
            const word& Sname
        ),
        (mesh, transportProperties, Sname)
    );


    // Selectors

    //- Return a reference to the selected relativePermeability model
    static autoPtr<relativePermeabilityModel> New
    (
        const fvMesh& mesh,
        const dictionary& transportProperties,
        const word& Sname
    );


    // Constructors

    //- Construct from components
    relativePermeabilityModel
    (
        const fvMesh& mesh,
        const dictionary& transportProperties,
        const word& Sname
    );


    //- Destructor
    virtual ~relativePermeabilityModel()
        {}


    // Member Functions

    //- Return the phase transport properties dictionary
    const dictionary& transportProperties() const
        {
            return transportProperties_;
        }

    //- Return relative permeabilities for phase a and b
    tmp<volScalarField> kra() const
        {
            return kra_;
        }

    tmp<volScalarField> krb() const
        {
            return krb_;
        }

    //- Return relative permeabilities for fractures
    tmp<volScalarField> krf() const
        {
            return krb_;
        }

    //- Return relative permeablities for matrix
    tmp<volScalarField> krm() const
        {
            return kra_;
        }

    //- Return relative permeability (one cell)
    scalar krb(const label& celli) const
        {
            return krb_[celli];
        }

    //- Return the derivatives of relative permeabilities
    tmp<volScalarField> dkradS() const
        {
            return dkradS_;
        }   
    tmp<volScalarField> dkrbdS() const
        {
            return dkrbdS_;
        }

    //- Return effective saturation
    tmp<volScalarField> Se() const
        {
            return Se_;
        }

    //- Set specific value of krb
    void setKrb(const label& celli, const scalar& krValue)
        {
            krb_[celli] = krValue;
        }

    //- Correct the relative permeabilities
    virtual void correct(const volScalarField& Sb) = 0;
    virtual void correctkra(const volScalarField& Sb) = 0;
    virtual void correctkrb(const volScalarField& Sb) = 0;
    virtual void correctkrb(const volScalarField& Sb, const label& celli) = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
