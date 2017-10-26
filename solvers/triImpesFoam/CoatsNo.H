// --- Coats (2003) IMPES STABILITY: CFL LIMIT
// - inertial part
{
    volScalarField dfw (
        "dfw",
        (dkrwdS*kra - dkradS*krw) /(Mmu*Foam::pow(kra,2)+2*kra*krw+1/Mmu*Foam::pow(krw,2))
    );

    dimensionedScalar smallRate("smallRate",dimVolume/dimTime, SMALL);

    dfw -= K*(rhoa-rhow)*fvc::surfaceSum(mag(mesh.Sf() & g))/fvc::surfaceSum(mag(phi)+smallRate)* (Foam::pow(kra,2)*dkrwdS/mua + Foam::pow(krw,2)*dkradS/muw)/(Mmu*Foam::pow(kra,2)+2*kra*krw+1/Mmu*Foam::pow(krw,2));

    scalarField CFLCoats((runTime.deltaT()/eps/mesh.V())*dfw*fvc::surfaceSum(mag(phi)));

    // - capillarity part
    if(activateCapillarity)
    {
        CFLCoats += (runTime.deltaT()/eps/mesh.V())*2*mag(pcModel->dpcawdS())*fvc::surfaceSum(Kf*mesh.magSf()/mag(mesh.delta()))*(kra*krw/(muw*kra+mua*krw));
    }

    Info<< "Coats CFL Number mean: " << gAverage(CFLCoats) << " max: " << gMax(CFLCoats) << endl;

    //CFLUse = gMax(CFLCoats);
    maxDeltaTFact = maxCo/(gMax(CFLCoats) + SMALL);
      
}
