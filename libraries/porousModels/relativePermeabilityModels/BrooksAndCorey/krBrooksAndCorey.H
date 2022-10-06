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
    Foam::relativePermeabilityModels::krBrooksAndCorey

Description
    Standard Brooks and Corey relative permeability model.

SourceFiles
    krBrooksAndCorey.C

\*---------------------------------------------------------------------------*/

#ifndef krBrooksAndCorey_H
#define krBrooksAndCorey_H

#include "relativePermeabilityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace relativePermeabilityModels
{

/*---------------------------------------------------------------------------*\
    Class krBrooksAndCorey Declaration
\*---------------------------------------------------------------------------*/

class krBrooksAndCorey
:
    public relativePermeabilityModel
{
    const dictionary& krBrooksAndCoreyCoeffs_;
    volScalarField n_;

    //- end points
    volScalarField kramax_;
    volScalarField krbmax_;
        
public:

    //- Runtime type information
    TypeName("BrooksAndCorey");

    // Constructors

    //- Construct from components
    krBrooksAndCorey
    (
        const fvMesh& mesh,
        const dictionary& transportProperties,
        const word& Sname,
        const word mediumName = ""
    );

    //- Destructor
    ~krBrooksAndCorey()
        {}

    // Member Functions

    //- Correct the relative permeabilities
    void correct(const volScalarField& Sb)
        {
            Se_= (Sb-Smin_)/(Smax_-Smin_);
            kra_ = kramax_* pow((scalar(1)-Se_),n_);
            krb_ = krbmax_ * pow(Se_,n_);
            dkradS_ = -kramax_*n_*pow((scalar(1)-Se_),n_-1)/(Smax_- Smin_);
            dkrbdS_ = krbmax_*n_*pow(Se_,n_-1)/(Smax_- Smin_);	
        }
    void correctkra(const volScalarField& Sb)
        {
            Se_= (Sb-Smin_)/(Smax_-Smin_);
            kra_ = kramax_* pow((scalar(1)-Se_),n_);
            dkradS_ = -kramax_*n_*pow((scalar(1)-Se_),n_-1)/(Smax_- Smin_);
        }
    void correctkrb(const volScalarField& Sb)
        {
            Se_= (Sb-Smin_)/(Smax_-Smin_);
            krb_ = krbmax_ * pow(Se_,n_);
            dkradS_ = -kramax_*n_*pow((scalar(1)-Se_),n_-1)/(Smax_- Smin_);
        }
    void correctkrb(const volScalarField& Sb, const label& celli)
        {
            scalar Se = (Sb[celli]-Smin_[celli])/(Smax_[celli]-Smin_[celli]);
            krb_[celli] = krbmax_[celli] * pow(Se,n_[celli]);
        }
    tmp<volScalarField> kr(const volScalarField& S)
        {
            volScalarField Se= (S-Smin_)/(Smax_-Smin_);
            return krbmax_ * pow(Se,n_);
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace relativePermeabilityModels

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
