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
#include "volFields.H"

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

    //- Field to write
    PtrList<const volScalarField> fieldsToWrite_;
    PtrList<const volScalarField> coeffFields_;
    PtrList<const surfaceScalarField> phiFields_;

    //- OFstream for CSV outputs
    PtrList<OFstream> CSVoutputFiles_;

    //- Internal function
    scalar computeInterpolationFactor();

public:

    static autoPtr<outputEventFile> New(Time&);

    // Constructors
    //- Construct from fileName
    outputEventFile
    (
        const word& fileName,
        const bool& isPresent,
        const bool& CSVoutput,
        Time&
    );

    // Selectors

    //- Destructor
    ~outputEventFile();

    // Access function
    bool isPresent() const{ return isPresent_;};

    // Member Functions
    scalar timeInterpolate(const scalar&, const scalar&);
    volScalarField timeInterpolate(const volScalarField&, bool = true);
    volVectorField timeInterpolate(const volVectorField&, bool = true);
    surfaceScalarField timeInterpolate(const surfaceScalarField&, bool = true);
    void checkControlDict() const;
    void addField(const volScalarField&, const volScalarField&, const surfaceScalarField&);
    void write();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
