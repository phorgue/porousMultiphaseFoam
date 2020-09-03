//- Reading MNT file (for seepage) if present
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

//- initialization of seepage option
bool seepageActive = transportProperties.lookupOrDefault("seepage", false);
if (seepageActive) Info << nl << "Seepage option is active" << endl;
else Info << nl << "no Seepage option" << endl;

if (seepageActive)
{
    if (potentialMNT.headerOk()) Info << nl << "Reading precomputed potentialMNT file in constant/" << endl;
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
        else FatalErrorIn("readFixedPoints.H") << nl << "no potentialMNT file neither MNT file while seepage is active " << endl;
    }

    //- checking that potentialMNT in superior to z0)
    if (min(potentialMNT.internalField() - z0.internalField()).value() <= 0)
    {
        Warning() << "potential MNT inferior to z0 in domain " << endl;
        forAll(potentialMNT.internalField(),celli)
        {
            potentialMNT.ref()[celli] = max(potentialMNT.internalField()[celli],z0.internalField()[celli]+0.01);
        }
    }
}

//- reading fixed potential list
List<Tuple2<point,scalar> > fixedPotentialList(transportProperties.lookupOrDefault
("fixedPotentialList",List<Tuple2<point,scalar> >())
);
labelList fixedPotentialIDList(fixedPotentialList.size());
scalarList fixedPotentialValueList(fixedPotentialList.size());
if (fixedPotentialList.size() > 0)
{
    forAll(fixedPotentialList,pointi)
    {
        fixedPotentialIDList[pointi] = mesh.findNearestCell(fixedPotentialList[pointi].first());
        fixedPotentialValueList[pointi] = fixedPotentialList[pointi].second();
    }

    //- Display information about fixed values
    Info << nl << "Fixed potential positions and values are " << nl << "{";
    forAll(fixedPotentialList,pointi)
    {
        Info << nl << "  " << fixedPotentialList[pointi].first() << " " << fixedPotentialValueList[pointi];
    }
    Info << nl << "}" << endl;
}

//- creating seepage and dryCell lists
labelList dryCellIDList(0);
labelList seepageIDList(0);
scalarList seepageValueList(0);
volScalarField cellFlux(fvc::div(phi*fvc::interpolate(hwater)));
