//- Reading DEM file name and potentialDEM if present
word DEMfileName(transportProperties.lookupOrDefault<word>("fileDEM",""));
volScalarField potentialDEM
(
    IOobject
    (
        "potentialDEM",
        runTime.constant(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    z0
);

//- Computing potentialDEM if necessary
if (potentialDEM.headerOk()) Info << nl << "Reading precomputed potentialDEM file in constant/" << endl;
else
{
    if (DEMfileName.size() > 0)
    {
        Info << nl << "Reading DEM file to compute potentialDEM...";
        DEMfile potentialDEMfile(DEMfileName);
        Info << "OK" << endl;

        //- Computing for potentialDEM for Seep term
        Info << "Interpolating value for potentialDEM...";
        forAll(mesh.C(),celli)
        {
            potentialDEM[celli] = potentialDEMfile.interpolate(mesh.C()[celli]);
        }
        Info << "OK" << endl;
        potentialDEM.write();
    }
    else Info << nl << "no potentialDEM field";
}

//- checking that potentialDEM in superior to z0)
if (gMin((potentialDEM.internalField() - z0.internalField())()) < 0)
{
    Warning() << "potential DEM inferior to z0 in domain => set to z0+hwaterMin" << endl;
    forAll(potentialDEM.internalField(),celli)
    {
        potentialDEM.ref()[celli] = max(potentialDEM.internalField()[celli],z0.internalField()[celli]+hwaterMin.value());
    }
}

//- initialization of seepage option
bool seepageActive = transportProperties.lookupOrDefault("seepage", false);
if (seepageActive)
{
    Info << nl << "Seepage option is active" << endl;
    if (!potentialDEM.headerOk() && DEMfileName.size() == 0)
    {
        FatalErrorIn("readFixedPoints.H") << nl << "no potentialDEM file neither DEM file while seepage is active " << abort(FatalError);
    }
}
else Info << nl << "no Seepage option" << endl;

//- reading fixed potential list
List<Tuple2<point,scalar> > fixedPotentialList(transportProperties.lookupOrDefault
("fixedPotentialList",List<Tuple2<point,scalar> >())
);
DynamicList<label> fixedPotentialIDList(0);
DynamicList<scalar> fixedPotentialValueList(0);
bool useDEMtoFixPotential(transportProperties.lookupOrDefault<bool>("useDEMtoFixPotential",false));
if (fixedPotentialList.size() > 0)
{
    label nFixedPoints = 0;
    //- find closest cell to fixed points
    forAll(fixedPotentialList,pointi)
    {
        label cellID =  mesh.findCell(fixedPotentialList[pointi].first());
        if (cellID > -1)
        {
            nFixedPoints += 1;
            fixedPotentialIDList.resize(nFixedPoints, cellID);
            fixedPotentialValueList.resize(nFixedPoints, fixedPotentialList[pointi].second());
        }
    }

    if (useDEMtoFixPotential)
    {
        Info << nl << "potentialDEM used for fixedPotentialList" << endl;
        forAll(fixedPotentialList,pointi)
        {
            fixedPotentialValueList[pointi] = potentialDEM[fixedPotentialIDList[pointi]];
        }
    }
    else
    {
        Info << nl << "user-defined values used for fixedPotentialList" << endl;
    }

    //- Display information about fixed values
    if (Pstream::nProcs() == 1)
    {
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
}

//- creating seepage and dryCell lists
labelList dryCellIDList(0);
labelList seepageIDList(0);
scalarList seepageValueList(0);
volScalarField cellFlux(fvc::div(phi*fvc::interpolate(hwater)));
