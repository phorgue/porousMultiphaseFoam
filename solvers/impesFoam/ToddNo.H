// --- Todd et al. (1972) Methods for increased accuracy in numerical reservoir simulator
{
    
    scalar count = 0.0;
    scalar somme = 0.0;
    scalarField Sum(Sb.size(),0);
    
    forAll( mesh.C() , celli) 
    {
        if(Sb[celli] > Sbmin.value() )
        {
            for(int i=0; i< mesh.cellCells(celli).size() ; i++)
            {
                if((Sb[celli] - Sb[ mesh.cellCells(celli)[i] ]) >= 0)
                {
                    Sum[celli] += (Sb[celli] - Sb[ mesh.cellCells(celli)[i] ]);}
                else
                {
                    Sum[celli] += ( Sb[ mesh.cellCells(celli)[i] ] - Sb [ celli ]);}
	 
                if(Sum[celli] < 0) { Info << "neg value !! " << Sb[ celli ] << ";" << Sb[ mesh.cellCells(celli)[i] ]<< abort(FatalError); }
	 
            }
            Sum[celli] /=  mesh.cellCells(celli).size() ;
            if(Sum[celli] > 0) { count++;}
        }
      
    }

    somme = gSum(Sum);
    somme /= max(count,1);


    if (runTime.controlDict().found("prefactor"))
    {
        scalar prefactor = runTime.controlDict().lookupOrDefault<scalar>("prefactor",-1);
        scalarField DS = Sb.oldTime().internalField() - Sb.internalField();
        scalar Todd_factor_kr = prefactor*somme/max(gMax(mag(DS)),SMALL);
        Todd_factor_kr = min(Todd_factor_kr,1e4);        
        Info << "Todd CFL kr = " << 1/Todd_factor_kr;
        maxDeltaTFact = Todd_factor_kr;
    }
    else
    {
        FatalErrorIn("in ToddNo.H")
            << "field 'prefactor' required in controlDict for Todd CFL condition : "
                << nl << "1D simulation : 1pt upwind => 0.5 , 2pt upwind => 0.25"
                << nl << "2D simulation : 1pt upwind => 0.25 , 2pt upwind => 0.125"
                << nl << "3D simulation : 1pt upwind => 0.1666 , 2pt upwind => 0.08333"       
                << abort(FatalError);
    }

    
    dimensionedScalar smallqP("smallRate",dimVolume/dimTime/p.dimensions(), VSMALL);
    if(activateCapillarity)
    {
        scalarField Tpc = mesh.V()/(mag(pcModel->dpcdS())*fvc::surfaceSum(Fbf*Maf*mesh.magSf()/mag(mesh.delta())
        + smallqP));
        scalar Todd_factor_pc = gMin(Tpc);

        Info << " , pc = " << gMin(Tpc) << endl;
        maxDeltaTFact = min(maxDeltaTFact,maxCo/Todd_factor_pc);
    }
    
}
