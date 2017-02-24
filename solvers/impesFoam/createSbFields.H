//- relative permeability (kr)
krModel->correct();
volScalarField kra = krModel->kra();
volScalarField krb = krModel->krb();
volScalarField dkradS = krModel->dkradS();
volScalarField dkrbdS = krModel->dkrbdS();


surfaceScalarField kraf ("kraf",fvc::interpolate(kra,"kra"));
surfaceScalarField krbf ("krbf",fvc::interpolate(krb,"krb"));
surfaceScalarField dkrafdS ("dkrafdS",fvc::interpolate(dkradS,"kra"));
surfaceScalarField dkrbfdS ("dkrbfdS",fvc::interpolate(dkrbdS,"krb"));

//- mobility and fractional flow 
surfaceScalarField Maf ("Maf",Kf*kraf/mua);
surfaceScalarField Laf ("Laf",rhoa*Kf*kraf/mua);	
surfaceScalarField Mbf ("Mbf",Kf*krbf/mub);
surfaceScalarField Lbf ("Lbf",rhob*Kf*krbf/mub);
surfaceScalarField Mf ("Mf",Maf+Mbf);
surfaceScalarField Lf ("Lf",Laf+Lbf);
surfaceScalarField Fbf ("Fbf",Mbf/Mf);
volScalarField Fb ("Fb",(krb/mub) / ( (kra/mua) + (krb/mub) ));

//- capillary pressure (pc)
myCapModel().correct();
surfaceVectorField gradpc("gradpc",fvc::interpolate(myCapModel().dpcdS()*fvc::grad(Sb),"pc"));

//- fluxes depending on saturation
surfaceScalarField phiG("phiG",(Lf * g) & mesh.Sf());
surfaceScalarField phiPc("phiPc",Mbf*gradpc & mesh.Sf());
