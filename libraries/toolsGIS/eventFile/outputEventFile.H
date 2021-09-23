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
    Foam::outputEventFile

Description
    Event file handler/reader for the porousMultiphaseFoam toolbox which
    contains 'm' dates used for time output instead of classical openfoam 
    write management. The class reads the file

    time1
    time2
    ...
    timen

SourceFiles
    outputEventFile.C

\*---------------------------------------------------------------------------*/

#ifndef outputEventFile_H
#define outputEventFile_H

#include "fvMesh.H"
#include "eventFile.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------* \
    Class fluidoutputEventFile Declaration
    \*---------------------------------------------------------------------------*/

class outputEventFile
:
    public eventFile
{
private:

    scalar computeInterpolationFactor(const Time&);

public:

    // Constructors    
    //- Copy constructor
    outputEventFile
    (
        const outputEventFile& fileToCopy
    );

    //- Construct from Istream
    outputEventFile
    (
        const word& fileName
    );
    
    // Selectors

    //- Destructor
    ~outputEventFile();

    // Member Functions
    scalar timeInterpolate(const scalar&, const scalar&, const Time&);
    volScalarField timeInterpolate(const volScalarField&, Time&, bool = true);
    volVectorField timeInterpolate(const volVectorField&, Time&, bool = true);
    surfaceScalarField timeInterpolate(const surfaceScalarField&, Time&, bool = true);
    void checkControlDict(const Time&) const;
    //- Access function
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
