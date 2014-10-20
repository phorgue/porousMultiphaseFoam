// --- Coats (2003) IMPES STABILITY: CFL LIMIT
// - inertial part


volTensorField dfw (
		    "dfw",
		    (dkrbdS*kra - dkradS*krb) /(Mmu*Foam::pow(kra,2)+2*kra*krb+1/Mmu*Foam::pow(krb,2))*tensor::I
			);

dimensionedScalar smallRate("smallRate",dimVolume/dimTime, SMALL);


// - gravity contribution
dfw -= (rhoa-rhob)*fvc::surfaceSum(Kf*mag(mesh.Sf() & g))/fvc::surfaceSum(mag(phi)+smallRate)* (Foam::pow(kra,2)*dkrbdS/mua + Foam::pow(krb,2)*dkradS/mub)/(Mmu*Foam::pow(kra,2)+2*kra*krb+1/Mmu*Foam::pow(krb,2));


volTensorField CFLCoats((runTime.deltaT()/eps/cellVolume)*dfw*fvc::surfaceSum(mag(phi)));


// - capillarity part of CFL
volScalarField dpc(myCapModel->dpcdS());
if(activateCapillarity)
  {
    CFLCoats += (runTime.deltaT()/eps/cellVolume)*2*mag(dpc)*fvc::surfaceSum(Kf*mesh.magSf()/mag(mesh.delta()))*(kra*krb/(mub*kra+mua*krb));  
  }


scalar TMP0( CFLCoats[0].component(0) );

forAll(mesh.C(),celli)
{
  scalar tmp0 =  CFLCoats[celli].component(0);
  
  for(int i=1; i<9 ; i++)
      {
	if(CFLCoats[celli].component(i)>tmp0) tmp0 = CFLCoats[celli].component(i);
      }
  if(TMP0<tmp0) TMP0 = tmp0;

}

Info<< "Coats CFL Number max: " << TMP0 << endl;
