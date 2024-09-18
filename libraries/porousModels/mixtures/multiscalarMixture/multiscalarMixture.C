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

#include "multiscalarMixture.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(multiscalarMixture, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multiscalarMixture::multiscalarMixture
(
    const dictionary& dict,
    const wordList& speciesNames,
    const fvMesh& mesh,
    const word& phaseName,
    List<sourceEventFile*>* sourceEventFileRegistry,
    const word& sourceEventDtFieldNameOverride,
    const dimensionSet& sourceTermDimFactor,
    const word& porousRegion
)
    :
    basicMultiComponentMixture(dict, speciesNames, mesh, phaseName),
    R_(speciesNames.size()),
    lambdas_(speciesNames.size()),
    dispersionModels_(speciesNames.size()),
    sourceEvents_(speciesNames.size()),
    sourceTerms_(speciesNames.size()),
    Kd_(speciesNames.size()),
    rs_(speciesNames.size()),
    epsTotal_(speciesNames.size())
{
    forAll(speciesNames, speciesi)
    {

        Info<< nl << "******** Specie " << speciesNames[speciesi] << porousRegion << " ********" << endl;

        dictionary specieDict(dict.optionalSubDict(speciesNames[speciesi]));
        specieDict = specieDict.optionalSubDict(porousRegion);
        dictionary porousTransport(specieDict.subDict("porousTransport"));

        dimensionedScalar Kd(porousTransport.getOrDefault("Kd", dimensionedScalar(dimensionSet(-1,3,0,0,0), 0)));
        Kd_[speciesi] = Kd.value();
        dimensionedScalar rs(porousTransport.getOrDefault("rs", dimensionedScalar(dimensionSet(1,-3,0,0,0), 0)));
        rs_[speciesi] = rs.value();
        dimensionedScalar epsTotal(porousTransport.getOrDefault("epsTotal", dimensionedScalar(dimless, 1)));
        epsTotal_[speciesi] = epsTotal.value();

        lambdas_[speciesi].dimensions().reset(dimensionSet(0,0,-1,0,0));
        lambdas_[speciesi] = porousTransport.getOrDefault("lambda", dimensionedScalar(dimensionSet(0,0,-1,0,0), 0));

        if (porousRegion != "")
        {
            Info << nl << "porousTransport " << porousRegion << " parameters" << nl << "{" << endl;
        }
        else
        {
            Info << nl << "porousTransport parameters" << nl << "{" << endl;
        }
        Info << "    Kd : " << Kd.value() << endl;
        Info << "    rs : " << rs.value() << endl;
        Info << "    epsTotal : " << epsTotal.value() << endl;
        Info << "    lambda : " << lambdas_[speciesi].name() << " : " << lambdas_[speciesi].value() << endl;
        Info << "}" << endl;

        R_.set
                (
                        speciesi,
                        new volScalarField
                                (
                                        IOobject
                                                (
                                                        Y(speciesi).name() + porousRegion + "_R",
                                                        mesh.time().timeName(),
                                                        mesh,
                                                        IOobject::NO_READ,
                                                        IOobject::NO_WRITE
                                                ),
                                        mesh,
                                        dimless,
                                        // calculatedFvPatchField<scalar>::typeName
                                        zeroGradientFvPatchField<scalar>::typeName
                                )
                );

        //- creation of dispersion model
        dispersionModels_.set(speciesi, dispersionModel::New("DeffModel", specieDict, mesh));

        //Handle source events

        const dimensionSet dimSourceTerm = Y(speciesi).dimensions()*sourceTermDimFactor/dimTime;

        //- read source event if present
        if(specieDict.found("sourceEventFileTracer"))
        {
            word sourceEventFileName = specieDict.getOrDefault<word>("sourceEventFileTracer","");
            sourceEvents_.set(speciesi, new sourceEventFile(sourceEventFileName, true));

            const word& dtFieldName =
                    sourceEventDtFieldNameOverride.empty() ? Y(speciesi).name() : sourceEventDtFieldNameOverride;

            sourceEvents_.last().setTimeScheme(dtFieldName, mesh);
            sourceEvents_.last().setFieldDimensions(dimSourceTerm);

            //- report found event to caller
            if(sourceEventFileRegistry)
            {
                sourceEventFileRegistry->append(&sourceEvents_.last());
            }
            else
            {
                FatalErrorIn("multiscalarMixtureI.H")
                        << "sourceEventFileTracer used with an incompatible solver"
                        << abort(FatalError);
            }

        }
        else//- otherwise, create a zero source term of the appropiate dimensions
        {
            sourceTerms_.set
                    (
                            speciesi,
                            new volScalarField
                                    (
                                            IOobject
                                                    (
                                                            "zeroSourceTerm",
                                                            mesh.time().timeName(),
                                                            mesh,
                                                            IOobject::NO_READ,
                                                            IOobject::NO_WRITE
                                                    ),
                                            mesh,
                                            dimensionedScalar("zero", dimSourceTerm, 0)
                                    )

                    );
        }

    }

}

// * * * * * * * * * * * * * * * * Members * * * * * * * * * * * * * * * * * //

void Foam::multiscalarMixture::correct
        (
                const volVectorField& U,
                const volScalarField& theta
        )
{
    forAll(Y(), speciesi)
    {
        dispersionModels_[speciesi].correct(Y(speciesi), U, theta);
        R_[speciesi].primitiveFieldRef() = 1 + (1-epsTotal_[speciesi]) * rs_[speciesi] * Kd_[speciesi] / theta;
        if(auto event = sourceEvents_.get(speciesi))
        {
            sourceTerms_.set(speciesi, event->dtValuesAsField());
        }
    }
}

void Foam::multiscalarMixture::correct
        (
                const volVectorField& U,
                const volScalarField& saturation,
                const volScalarField& eps
        )
{

    forAll(Y(), speciesi)
    {
        dispersionModels_[speciesi].correct(Y(speciesi), U, saturation, eps);
        R_[speciesi].primitiveFieldRef() = 1 + (1-epsTotal_[speciesi]) * rs_[speciesi] * Kd_[speciesi] / (eps*saturation);
        if(auto event = sourceEvents_.get(speciesi))
        {
            sourceTerms_.set(speciesi, event->dtValuesAsField());
        }
    }
}

bool Foam::multiscalarMixture::initRetardCoef(const volScalarField& eps)
{
    Info << nl << "Computing initial retard coefficient: ";
    forAll(Y(), speciesi)
        {
            if (rs_[speciesi].value() != 0 &&  Kd_[speciesi].value() != 0)
            {
                if (gMax(eps) > 1)
                {
                    FatalErrorIn("multiscalarMixture.C") <<
                        "Field " << eps.name() << " seems to be uninitialized. You should specify value in transportProperties or field in constant/" << abort(FatalError);
                }
                R_[speciesi].primitiveFieldRef() = 1 + (1-epsTotal_[speciesi]) * rs_[speciesi] * Kd_[speciesi] / eps;
                Info << "Min(R) = " << gMin(R_[speciesi]) << " Max(R) = " << gMax(R_[speciesi]) << endl;
            }
        }
    return false;

}

// ************************************************************************* //
