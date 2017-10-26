//- relative permeability (kr)
krModel->correct();
volScalarField kra = krModel->kra();
volScalarField kro = krModel->kro();
volScalarField krw = krModel->krw();
volScalarField dkradS = krModel->dkradS();
volScalarField dkrwdS = krModel->dkrwdS();

surfaceScalarField kraf ("kraf",fvc::interpolate(kra,"kra"));
surfaceScalarField krof ("krof",fvc::interpolate(kro,"kro"));
surfaceScalarField krwf ("krwf",fvc::interpolate(krw,"krw"));

//- mobility and fractional flow 
surfaceScalarField Maf ("Maf",Kf*kraf/mua);
surfaceScalarField Laf ("Laf",rhoa*Kf*kraf/mua);
surfaceScalarField Mof ("Mof",Kf*krof/muo);
surfaceScalarField Lof ("Lof",rhoo*Kf*krof/muo);
surfaceScalarField Mwf ("Mwf",Kf*krwf/muw);
surfaceScalarField Lwf ("Lwf",rhow*Kf*krwf/muw);
surfaceScalarField Mf ("Mf",Maf+Mwf+Mof);
surfaceScalarField Lf ("Lf",Laf+Lwf+Lof);
// surfaceScalarField Fbf ("Fbf",Mbf/Mf);
volScalarField Fw ("Fb",(krw/muw) / ( (kra/mua) + (krw/muw) + (kro/muo) ));
volScalarField Fo ("Fb",(kro/muo) / ( (kra/mua) + (krw/muw) + (kro/muo) ));

//- capillary pressure (pc)
pcModel().correct();
//surfaceVectorField gradpc_o("gradpc_o",fvc::interpolate(pcModel().dpcaodS()*fvc::grad(So),"pc"));
//surfaceVectorField gradpc_w("gradpc_w", fvc::interpolate( (pcModel().dpcaodS()+pcModel().dpcowdS() )*fvc::grad(Sw),"pc"));

//Info << nl << (pcModel().dpcawdS() ) << endl;

//- fluxes depending on saturation
surfaceScalarField phiG("phiG",(Lf * g) & mesh.Sf());
surfaceScalarField phiPc_o("phiPc_o",Mof * fvc::interpolate(pcModel().dpcaodS(),"pc")* fvc::snGrad(So) * mesh.magSf());
//surfaceScalarField phiPc_w("phiPc_w",Mwf * fvc::interpolate(pcModel().dpcawdS()+pcModel().dpcowdS(),"pc")* fvc::snGrad(Sw) * mesh.magSf());
surfaceScalarField phiPc_w("phiPc_w",Mwf * fvc::interpolate(pcModel().dpcaodS()+pcModel().dpcowdS(),"pc")* fvc::snGrad(Sw) * mesh.magSf());

surfaceScalarField phiPc("phiPc",phiPc_o+phiPc_w);

// Info << nl << "VALIDATION"
// << nl <<  fvc::interpolate(pcModel().dpcaodS()+pcModel().dpcowdS(),"pc")
// << nl <<  fvc::interpolate(pcModel().dpcawdS(),"pc") << endl;
