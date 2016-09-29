// --- Coats (2003) IMPES STABILITY: CFL LIMIT
// - inertial part
{
volScalarField dfw (
		    "dfw",
		    (dkrbdS*kra - dkradS*krb) /(Mmu*Foam::pow(kra,2)+2*kra*krb+1/Mmu*Foam::pow(krb,2))
			);

dimensionedScalar smallRate("smallRate",dimVolume/dimTime, SMALL);

dfw -= K*(rhoa-rhob)*fvc::surfaceSum(mag(mesh.Sf() & g))/fvc::surfaceSum(mag(phi)+smallRate)* (Foam::pow(kra,2)*dkrbdS/mua + Foam::pow(krb,2)*dkradS/mub)/(Mmu*Foam::pow(kra,2)+2*kra*krb+1/Mmu*Foam::pow(krb,2));

volScalarField CFLCoats((runTime.deltaT()/eps/cellVolume)*dfw*fvc::surfaceSum(mag(phi)));

// - capillarity part
if(activateCapillarity)
  {
    CFLCoats += (runTime.deltaT()/eps/cellVolume)*2*mag(myCapModel->dpcdS())*fvc::surfaceSum(Kf*mesh.magSf()/mag(mesh.delta()))*(kra*krb/(mub*kra+mua*krb));  
  }

Info<< "Coats CFL Number mean: " << gAverage(CFLCoats) << " max: " << gMax(CFLCoats) << endl;

      CFLUse = gMax(CFLCoats.internalField());
      maxDeltaTFact = maxCo/(CFLUse + SMALL);
      
}
