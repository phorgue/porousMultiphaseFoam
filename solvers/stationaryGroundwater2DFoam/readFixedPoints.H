//- Read the global MNT mode : fixedMNT, fixed or seepage
word forcingPotential=transportProperties.lookupOrDefault<word>("forcingPotential","none");
Info << nl << "Forcing potential mode is : " << forcingPotential << endl;

//- Reading MNT file and computing potentialMNT field if required
volScalarField potentialMNT(
    IOobject
    (
        "potentialMNT",
        runTime.constant(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("default_potential",dimLength,0),
    zeroGradientFvPatchField<scalar>::typeName
);
if ((forcingPotential == "fixedMNT") || (forcingPotential == "seepage"))
{
    if (potentialMNT.headerOk())
    {
        Info << nl << "Reading precomputed potentialMNT file in constant/" << endl;
    }
    else
    {
        word MNTfileName(transportProperties.lookupOrDefault<word>("fileMNT",""));
        if (MNTfileName.size() > 0)
        {
            Info << nl << "Reading MNT file to compute potentialMNT...";
            MNTfile potentialMNTfile(MNTfileName);
            Info << "OK" << endl;

            //- Computing for potentialMNT for Seep term
            Info << "Interpolating value for potentialMNT...";
            forAll(mesh.C(),celli)
            {
                potentialMNT[celli] = potentialMNTfile.interpolate(mesh.C()[celli]);
            }
            Info << "OK" << endl;
            potentialMNT.write();
        }
        else
        {
            FatalErrorIn("readFixedPoints.H") << nl << "no potentialMNT file neither MNT file for forcing potential mode : "
                << forcingPotential << endl;
        }
    }

    //- checking that potentialMNT in superior to z0)
    if (min(potentialMNT.internalField() - z0.internalField()).value() <= 0)
    {
        Info << "*** WARNING *** : potential MNT inferior to z0 in domain " << endl;
        forAll(potentialMNT.internalField(),celli)
        {
            potentialMNT.ref()[celli] = max(potentialMNT.internalField()[celli],z0.internalField()[celli]+0.01);
        }
    }
}
else if ((forcingPotential != "fixedValue") && (forcingPotential != "none"))
{
    FatalErrorIn("readFixedPoint.H") << "forcingPotential word unknown : " << forcingPotential << abort(FatalError);
}

List<Tuple2<point,scalar> > fixedPotentialList(transportProperties.lookupOrDefault
("fixedPotentialList",List<Tuple2<point,scalar> >())
);
labelList fixedPotentialIDList(fixedPotentialList.size());
scalarList fixedPotentialValueList(fixedPotentialList.size());
volScalarField seepTerm = 0*potentialMNT/runTime.deltaT();

if (forcingPotential.substr(0,5) == "fixed")
{
    //- find nearest cell
    forAll(fixedPotentialList,pointi)
    {
        fixedPotentialIDList[pointi] = mesh.findNearestCell(fixedPotentialList[pointi].first());
    }

    if (forcingPotential == "fixedValue")
    {
        forAll(fixedPotentialList,pointi)
        {
            fixedPotentialValueList[pointi] = fixedPotentialList[pointi].second();
        }
    }
    else if (forcingPotential == "fixedMNT")
    {
        forAll(fixedPotentialList,pointi)
        {
            fixedPotentialValueList[pointi] =  potentialMNT[fixedPotentialIDList[pointi]];
        }
    }

    //- Display information about fixed values
    Info << nl << "Fixed potential positions and values are " << nl << "{";
    forAll(fixedPotentialList,pointi)
    {
        Info << nl << "  " << fixedPotentialList[pointi].first() << " " << fixedPotentialValueList[pointi];
    }
    Info << nl << "}" << endl;
}

volScalarField cellFlux(fvc::div(phi*fvc::interpolate(hwater)));
