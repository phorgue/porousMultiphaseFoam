/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.
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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<typename EpsType>
Foam::multiscalarMixture::multiscalarMixture
(
    const dictionary& dict,
    const wordList& speciesNames,
    const fvMesh& mesh,
    const word& phaseName,
    const EpsType& eps,
    List<sourceEventFile*>* sourceEventFileRegistry,
    const word& sourceEventDtFieldNameOverride,
    const dimensionSet& dimSourceTerm
)
:
    basicMultiComponentMixture(dict, speciesNames, mesh, phaseName),
    R_(speciesNames.size()),
    lambdas_(speciesNames.size()),
    dispersionModels_(speciesNames.size()),
    sourceEvents_(speciesNames.size()),
    zeroSourceTerm_
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
{
    forAll(speciesNames, speciesi)
    {

        Info<< "Species " << speciesNames[speciesi] << endl;

        dictionary speciesDict(dict.optionalSubDict(speciesNames[speciesi]));
        dictionary porousTransport(speciesDict.subDict("porousTransport"));

        dimensionedScalar Kd(porousTransport.lookupOrDefault("Kd",dimensionedScalar("Kd_default",dimensionSet(-1,3,0,0,0),0.)));
        dimensionedScalar rs(porousTransport.lookupOrDefault("rs",dimensionedScalar("rs_default",dimensionSet(1,-3,0,0,0),0.)));
        dimensionedScalar epsTotal(porousTransport.lookupOrDefault("epsTotal",dimensionedScalar("epsTotal_default",dimless,1.)));

        lambdas_[speciesi].dimensions().reset(dimensionSet(0,0,-1,0,0));
        lambdas_[speciesi] = porousTransport.lookupOrDefault("lambda",dimensionedScalar("lambda_default",dimensionSet(0,0,-1,0,0),0.));

        Info << nl << "porousTransport parameters" << nl << "{" << endl;
        Info << "    " << Kd.name() << " : " << Kd.value() << endl;
        Info << "    " << rs.name() << " : " << rs.value() << endl;
        Info << "    " << epsTotal.name() << " : " << epsTotal.value() << endl;
        Info << "    " << lambdas_[speciesi].name() << " : " << lambdas_[speciesi].value() << endl;
        Info << "}" << endl;

        Info << nl << "Computing retard coefficient: ";

        R_.set
        (
            speciesi,
            new volScalarField
            (
                IOobject
                (
                    Y(speciesi).name() + "_R",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimless
            )
        );

        R_[speciesi] =  1 + (1-epsTotal) * rs * Kd / eps;
        Info << "Min(R) = " << gMin(R_[speciesi]) << " Max(R) = " << gMax(R_[speciesi]) << endl;

        //- creation of dispersion model
        dispersionModels_.set(speciesi, dispersionModel::New("DeffModel", speciesDict, mesh));


        //- read source event if present
        if(speciesDict.found("eventFileTracerSource"))
        {
            word sourceEventFileName = speciesDict.lookupOrDefault<word>("eventFileTracerSource","");
            sourceEvents_.set(speciesi, new sourceEventFile(sourceEventFileName));

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
                FatalErrorIn("multiscalarMixture.C")
                    << "eventFileTracerSource used with an incompatible solver"
                    << abort(FatalError);
            }

        }

    } 

}


// * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * //


inline const Foam::volScalarField& Foam::multiscalarMixture::R(const label speciesi) const
{
    return R_[speciesi];
}


inline const Foam::dimensionedScalar& Foam::multiscalarMixture::lambda(const label speciesi) const
{
    return lambdas_[speciesi];
}


template<typename ThetaType>
Foam::tmp<Foam::volTensorField> Foam::multiscalarMixture::Deff
(
    const label speciesi,
    const volVectorField& U,
    const ThetaType& theta
) const
{
    auto& dispersionModel = dispersionModels_[speciesi];

    dispersionModel.correct(Y(speciesi), U, theta);

    return dispersionModel.Deff();
}


template<typename SaturationType, typename EpsType>
Foam::tmp<Foam::volTensorField> Foam::multiscalarMixture::Deff
(
    const label speciesi,
    const volVectorField& U,
    const SaturationType& saturation,
    const EpsType& eps
) const
{
    auto& dispersionModel = dispersionModels_[speciesi];

    dispersionModel.correct(Y(speciesi), U, saturation, eps);

    return dispersionModel.Deff();
}


inline Foam::tmp<Foam::volScalarField> Foam::multiscalarMixture::sourceTerm(const label speciesi) const
{
    if (sourceEvents_(speciesi))
    {
        return sourceEvents_[speciesi].dtValuesAsField();
    }
    else return zeroSourceTerm_;
}



// ************************************************************************* //