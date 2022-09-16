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
    Foam::relativePermeabilityModels::krIppisch

Description
    Ippisch Vogel Bastian relative permeability model (AWR 2006)

SourceFiles
    krIppisch.C

\*---------------------------------------------------------------------------*/

#ifndef krIppisch_H
#define krIppisch_H

#include "relativePermeabilityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace relativePermeabilityModels
{

/*---------------------------------------------------------------------------* \
    Class krIppisch Declaration
    \*---------------------------------------------------------------------------*/

class krIppisch
:
    public relativePermeabilityModel
{
    //- Private data
    volScalarField Smin_;
    volScalarField Smax_;

    //- Van Genuchten coefficients
    dictionary krIppischCoeffs_;
    volScalarField m_,n_,alpha_,tau_,he_;

    //- cut-off point saturation (computed with entry pressure)
    volScalarField Sc_;

public:

    //- Runtime type information
    TypeName("Ippisch");

    // Constructors

    //- Construct from components
    krIppisch
    (
        const fvMesh& mesh,
        const dictionary& transportProperties,
        const word& Sname
    );

    //- Destructor
    ~krIppisch()
        {}

    // Member Functions

    //- Correct the relative permeabilities
     void correct(const volScalarField& Sb)
        {
            FatalErrorIn("krIppisch.H") << " Ippisch model can be used only with groundwater solvers" << abort(FatalError);
        }
    void correctkra(const volScalarField& Sb)
        {
            FatalErrorIn("krIppisch.H") << " Ippisch model can be used only with groundwater solvers" << abort(FatalError);
        }
    void correctkrb(const volScalarField& Sb)
        {
            Se_= (Sb-Smin_)/(Smax_-Smin_);
            krb_ = pow(Se_,tau_) * pow(
                (1-pow(1-pow(Sc_*Se_,1/m_),m_)) /
                (1-pow(1-pow(Sc_,1/m_),m_))
                ,2);
        }
    void correctkrb(const volScalarField& Sb, const label& celli)
        {
            scalar Se = (Sb[celli]-Smin_[celli])/(Smax_[celli]-Smin_[celli]);
            krb_[celli] = pow(Se_[celli],tau_[celli]) * pow(
                (1-pow(1-pow(Sc_[celli]*Se,1/m_[celli]),m_[celli])) /
                (1-pow(1-pow(Sc_[celli],1/m_[celli]),m_[celli]))
                ,2);
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace relativePermeabilityModels

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
