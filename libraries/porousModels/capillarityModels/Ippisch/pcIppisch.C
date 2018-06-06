/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "pcIppisch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace capillarityModels
{
defineTypeNameAndDebug(pcIppisch, 0);

addToRunTimeSelectionTable
(
    capillarityModel,
    pcIppisch,
    dictionary
);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::capillarityModels::pcIppisch::pcIppisch
(
    const word& name,
    const dictionary& transportProperties,
    const volScalarField& Sb
)
    :
    capillarityModel(name, transportProperties,Sb),
    pcIppischCoeffs_(transportProperties.subDict(typeName + "Coeffs")),
    Smin_
    (
        IOobject
        (
            Sb_.name()+"min",
            Sb_.time().timeName(),
            Sb_.db(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        Sb.mesh(),
        transportProperties.lookupOrDefault(Sb_.name()+"min",dimensionedScalar(Sb_.name()+"min",dimless,0))
    ),
    Smax_
    (
        IOobject
        (
            Sb_.name()+"max",
            Sb_.time().timeName(),
            Sb_.db(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        Sb.mesh(),
        transportProperties.lookupOrDefault(Sb_.name()+"max",dimensionedScalar(Sb_.name()+"min",dimless,0))
    ),
    m_
    (
        IOobject
        (
            "m",
            Sb_.time().timeName(),
            Sb_.db(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        Sb.mesh(),
        dimensionedScalar("m",dimless,pcIppischCoeffs_.lookupOrDefault<scalar>("m",0))
    ),
    n_(1/(1-m_)),
    alpha_
    (
        IOobject
        (
            "alpha",
            Sb_.time().timeName(),
            Sb_.db(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        Sb.mesh(),
        dimensionedScalar("alpha",dimless,pcIppischCoeffs_.lookupOrDefault<scalar>("alpha",GREAT))
    ),
    tau_
    (
        IOobject
        (
            "tau",
            Sb_.time().timeName(),
            Sb_.db(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        Sb.mesh(),
        dimensionedScalar("tau",dimless,pcIppischCoeffs_.lookupOrDefault<scalar>("tau",1.))
    ),
    he_
    (
        IOobject
        (
            "he",
            Sb_.time().timeName(),
            Sb_.db(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        Sb.mesh(),
        dimensionedScalar("he",dimless,pcIppischCoeffs_.lookupOrDefault<scalar>("he",1.))
    ),
    Se_((Sb_-Smin_)/(Smax_-Smin_)),
    Sc_(pow(1+pow(alpha_*he_,n_),-m_))
{
    if (gMin(m_) == 0) FatalErrorIn("Foam::capillarityModels::pcIppisch::pcIppisch") << "m = 0 in pcIppisch" << abort(FatalError);
    Info << "Ippisch parameters for capillary pressure model" << nl << "{" << endl;
    Info << "    m ";
    if (m_.headerOk()) { Info << "read file" << endl;}
    else {Info << average(m_).value() << endl;}
    Info <<  "    n ";
    if (n_.headerOk()) { Info << "read file" << endl;}
    else {Info << average(n_).value() << endl;}
        Info <<  "    alpha ";
    if (alpha_.headerOk()) { Info << "read file" << endl;}
    else {Info << average(alpha_).value() << endl;}
        Info <<  "    tau ";
    if (tau_.headerOk()) { Info << "read file" << endl;}
    else {Info << average(tau_).value() << endl;}
    Info <<  "    he ";
    if (he_.headerOk()) { Info << "read file" << endl;}
    else {Info << average(he_).value() << endl;}
    Info << "} \n" << endl;     
}

// ************************************************************************* //

