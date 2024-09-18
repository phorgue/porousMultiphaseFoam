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

#include "zeroGradientFvPatchField.H"

// * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * //

inline const Foam::volScalarField& Foam::multiscalarMixture::R(const label speciesi) const
{
    return R_[speciesi];
}

inline const Foam::dimensionedScalar& Foam::multiscalarMixture::lambda(const label speciesi) const
{
    return lambdas_[speciesi];
}

inline const Foam::dimensionedScalar& Foam::multiscalarMixture::Dm(const label speciesi) const
{
    return dispersionModels_[speciesi].Dm();
}

inline const Foam::volTensorField& Foam::multiscalarMixture::Deff(const label speciesi) const
{
    return dispersionModels_[speciesi].Deff();
}

inline const Foam::volScalarField& Foam::multiscalarMixture::sourceTerm(const label speciesi) const
{
    return sourceTerms_[speciesi];
}

inline const Foam::PtrList<Foam::sourceEventFile>& Foam::multiscalarMixture::sourceEvents() const
{
    return sourceEvents_;
}

// ************************************************************************* //
