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
    Foam::dispersionModels::dualAlphaDispersion

Description
    Dispersion model based on alphaL and alphaT coefficients (Bear, 1979)

      Deff= Dm/tau + (alphaT * mag(U) + (alphaL - alphaT) * (U*U)/mag(U) )/theta

    where theta is water content :

      theta = Saturation * porosity

    Two sets of parameters are used for saturated (alphaL_sat/alphaT_sat) and
    unsaturated areas (alphaL_unsat/alphaT_unsat). Transition between values
    are computed using theta thresholds (thetaThreshold_unsat and
    thetaThreshold_sat)

SourceFiles
    dualAlphaDispersion.C

\*---------------------------------------------------------------------------*/

#ifndef dualAlphaDispersion_H
#define dualAlphaDispersion_H

#include "dispersionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace dispersionModels
{

/*---------------------------------------------------------------------------*\
                  Class dualAlphaDispersion Declaration
\*---------------------------------------------------------------------------*/

class dualAlphaDispersion
:
    public dispersionModel
{
    //- Private data

    //- Dictionary for model coefficients
    dictionary dualAlphaDispersionCoeffs_;

    //- model coefficients
    volScalarField tau_,alphaLsat_,alphaTsat_,alphaLunsat_,alphaTunsat_;
    volScalarField coef_, thetaThresholdSat_, thetaThresholdUnsat_;

public:

    //- Runtime type information
    TypeName("dualAlphaDispersion");

    // Constructors

    //- Construct from components
    dualAlphaDispersion
    (
        const word& name,
        const dictionary& transportProperties,
        const fvMesh& mesh
    );

    //- Destructor
    ~dualAlphaDispersion()
        {}

    // Member Functions

    //- Correct the dispersion giving water content or porosity field
    void correct(const volScalarField& C, const volVectorField& U, const volScalarField& theta)
        {
            dimensionedScalar USMALL("",dimLength/dimTime,1e-36);
            coef_ = min(1.0, max(0.0, (theta - thetaThresholdUnsat_) / (thetaThresholdSat_ - thetaThresholdUnsat_)));
            Deff_ = (Dm_/tau_ +  ((1-coef_) * alphaTunsat_ + coef_*alphaTsat_) * mag(U) / theta )*tensor::I
                + ( (1-coef_) * alphaLunsat_ + coef_*alphaLsat_ - (1-coef_) * alphaTunsat_ - coef_*alphaTsat_ ) * (U*U)/(theta*mag(U)+USMALL);

            forAll(Deff_.boundaryField(),patchi)
            {
                //- Dispersion equal to zero on patch where flux is fixed
                if (C.boundaryField().types()[patchi] == "eventFlux" ||
                    C.boundaryField().types()[patchi] == "inletOutlet"
                )
                {
                    Deff_.boundaryFieldRef()[patchi] = tensor(zero());
                }
            }
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace dispersionModels

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
