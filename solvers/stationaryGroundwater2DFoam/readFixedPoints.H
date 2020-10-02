//- Reading MNT file name and potentialMNT if present
word MNTfileName(transportProperties.lookupOrDefault<word>("fileMNT",""));
volScalarField potentialMNT
(
    IOobject
    (
        "potentialMNT",
        runTime.constant(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    z0
);

//- Computing potentialMNT if necessary
if (potentialMNT.headerOk()) Info << nl << "Reading precomputed potentialMNT file in constant/" << endl;
else
{
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
    else Info << nl << "no potentialMNT field";
}

//- checking that potentialMNT in superior to z0)
if (min(potentialMNT.internalField() - z0.internalField()).value() < 0)
{
    Warning() << "potential MNT inferior to z0 in domain => set to z0+hwaterMin" << endl;
    forAll(potentialMNT.internalField(),celli)
    {
        potentialMNT.ref()[celli] = max(potentialMNT.internalField()[celli],z0.internalField()[celli]+hwaterMin.value());
    }
}

//- initialization of seepage option
bool seepageActive = transportProperties.lookupOrDefault("seepage", false);
if (seepageActive)
{
    Info << nl << "Seepage option is active" << endl;
    if (!potentialMNT.headerOk() && MNTfileName.size() == 0)
    {
        FatalErrorIn("readFixedPoints.H") << nl << "no potentialMNT file neither MNT file while seepage is active " << abort(FatalError);
    }
}
else Info << nl << "no Seepage option" << endl;

//- reading fixed potential list
List<Tuple2<point,scalar> > fixedPotentialList(transportProperties.lookupOrDefault
("fixedPotentialList",List<Tuple2<point,scalar> >())
);
labelList fixedPotentialIDList(fixedPotentialList.size());
scalarList fixedPotentialValueList(fixedPotentialList.size());
bool useMNTtoFixPotential(transportProperties.lookupOrDefault<bool>("useMNTtoFixPotential",false));
if (fixedPotentialList.size() > 0)
{
    //- find closest cell to fixed points
    forAll(fixedPotentialList,pointi)
    {
        fixedPotentialIDList[pointi] = mesh.findNearestCell(fixedPotentialList[pointi].first());
    }

    if (useMNTtoFixPotential)
    {
        Info << nl << "potentialMNT used for fixedPotentialList" << endl;
        forAll(fixedPotentialList,pointi)
        {
            fixedPotentialValueList[pointi] = potentialMNT[fixedPotentialIDList[pointi]];
        }
    }
    else
    {
        Info << nl << "user-defined values used for fixedPotentialList" << endl;
        forAll(fixedPotentialList,pointi)
        {
            fixedPotentialValueList[pointi] = fixedPotentialList[pointi].second();
        }
    }

    //- Display information about fixed values
    Info << nl << "Fixed potential positions and values are " << nl << "{";
    forAll(fixedPotentialList,pointi)
    {
        scalar distance_to_centre = Foam::sqrt(pow(mesh.C()[fixedPotentialIDList[pointi]].x()-fixedPotentialList[pointi].first().x(),2)
        +pow(mesh.C()[fixedPotentialIDList[pointi]].y()-fixedPotentialList[pointi].first().y(),2));
        Info << nl << "  " << fixedPotentialList[pointi].first() <<  " : value " << fixedPotentialValueList[pointi]
            << "  (cellID = " << fixedPotentialIDList[pointi] << ", distance with cell-center = " << distance_to_centre << ")";
    }
    Info << nl << "}" << endl;
}

//- creating seepage and dryCell lists
labelList dryCellIDList(0);
labelList seepageIDList(0);
scalarList seepageValueList(0);
volScalarField cellFlux(fvc::div(phi*fvc::interpolate(hwater)));
