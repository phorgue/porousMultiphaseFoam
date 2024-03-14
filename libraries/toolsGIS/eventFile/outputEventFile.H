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

#include "eventFile.H"
#include "outputField.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------* \
    Class outputEventFile Declaration
\*---------------------------------------------------------------------------*/

class outputEventFile
:
    public eventFile
{
private:

    //- If output event file is present
    bool isPresent_;

    //- If CSV file output is active
    bool CSVoutput_;

    //- Link to Time
    Time& runTime_;

    //- Time interpolation factor
    scalar ifactor_;

    //- z-scale (for 2D solvers)
    scalar zscale_;

    //- One field for
    volScalarField one_;

    //- Field to write
    PtrList<outputField> outputFields_;

    //- Internal function
    void updateInterpolationFactor();

public:

    static autoPtr<outputEventFile> New(Time&, const fvMesh&, const scalar zscale = 1);

    // Constructors
    //- Construct from fileName
    outputEventFile
    (
        const word& fileName,
        const bool& isPresent,
        const bool& CSVoutput,
        Time&,
        const fvMesh&,
        const scalar zscale
    );

    // Selectors

    //- Destructor
    ~outputEventFile();

    // Access function
    bool isPresent() const{ return isPresent_;};
    bool CSVoutput() const{ return CSVoutput_;};
    scalar ifactor() const{ return ifactor_;};

    // Member functions
    void checkControlDict() const;
    void addField(
            const volScalarField& field,
            const surfaceScalarField& phi,
            const volScalarField& coef1,
            const volScalarField& coef2,
            const volScalarField& coef3,
            word fileName = "/dev/null",
            bool saturation = false
    );

    void addField(
            const volScalarField& field,
            const surfaceScalarField& phi,
            const volScalarField& coef1,
            const volScalarField& coef2,
            word fileName = "/dev/null",
            bool saturation = false
    );

    void addField(
            const volScalarField& field,
            const surfaceScalarField& phi,
            const volScalarField& coef1,
            word fileName = "/dev/null",
            bool saturation = false
    );

    void addField(
            const volScalarField& field,
            const surfaceScalarField& phi,
            word fileName = "/dev/null",
            bool saturation = false
    );
    void write();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
