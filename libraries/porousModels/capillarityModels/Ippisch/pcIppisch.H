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
    Foam::capillarityModels::Ippisch

Description
    Ippisch-Vogel-Bastian capillarity pressure model (AWR 2006)

SourceFiles
    pcIppisch.C

\*---------------------------------------------------------------------------*/

#ifndef pcIppisch_H
#define pcIppisch_H

#include "capillarityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace capillarityModels
{

/*---------------------------------------------------------------------------*\
                    Class pcIppisch Declaration
\*---------------------------------------------------------------------------*/

class pcIppisch
:
    public capillarityModel
{
    //- Ippisch coefficients
    volScalarField m_,n_,alpha_,tau_,he_;

    //- cut-off point saturation (computed with entry pressure)
    volScalarField Sc_;
    
public:

    //- Runtime type information
    TypeName("Ippisch");

    // Constructors

    //- Construct from components
    pcIppisch
    (
        const fvMesh& mesh,
        const dictionary& transportProperties,
        const word& Sname,
        const word porousRegion = ""
    );

    //- Destructor
    ~pcIppisch()
        {}

    // Member Functions

    //- Correct the capillarity pressure
    void correct(const volScalarField& Sb)
        {
            FatalErrorIn("void correct()")
                << "Not implemented in pcIppish.H " << abort(FatalError);
        }

    //- Compute capillarity capacity from head pressure and return moisture (Richards' model)
    tmp<volScalarField> correctAndSb(volScalarField& h)
        {
            //- Compute Saturation from head pressure
            h.dimensions().reset(dimless);
            Se_ = neg(h+he_) *(1/Sc_)*pow(1.0+pow(alpha_*mag(h),n_),-m_) + pos(h+he_+VSMALL) ;
            h.dimensions().reset(dimLength);

            //- Update capillarity capacity
            Ch_.dimensions().reset(dimless);
            Ch_ = alpha_*m_*(Smax_-Smin_)/(1.0-m_)*pow(Se_,1.0/m_)*pow(1.0-pow(Se_,1.0/m_),m_) ;
            Ch_.dimensions().reset(dimless/dimLength);

            //- Compute saturation
            return Smin_+(Smax_-Smin_)*Se_;
        }

    //- Compute capillarity capacity from head pressure and return moisture (Richards' model)
    scalar correctAndSb(volScalarField& h, const label& celli)
        {
            //- Compute saturation from head pressure
            scalar Se = neg(h[celli]+he_[celli]) *(1.0/Sc_[celli])*pow(1.0+pow(alpha_[celli]*mag(h[celli]),n_[celli]),-m_[celli]) + pos(h[celli]+he_[celli]+VSMALL) ;

            //- Update capillarity capacity
            Ch_[celli] = alpha_[celli]*m_[celli]*(Smax_[celli]-Smin_[celli])/(1.0-m_[celli])*pow(Se,1.0/m_[celli])*pow(1.0-pow(Se,1.0/m_[celli]),m_[celli]) ;

            //- Update saturation
            return Smin_[celli]+(Smax_[celli]-Smin_[celli])*Se;
        }

    //- Compute saturation from given head pressure
    tmp<volScalarField> S(volScalarField& h)
        {
            //- Compute Saturation from head pressure
            h.dimensions().reset(dimless);
            tmp<volScalarField> Se = neg(h+he_) *(1/Sc_)*pow(1.0+pow(alpha_*mag(h),n_),-m_) + pos(h+he_+VSMALL) ;
            h.dimensions().reset(dimLength);
            //- Return saturation
            return Smin_+(Smax_-Smin_)*Se;
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace capillarityModels

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
