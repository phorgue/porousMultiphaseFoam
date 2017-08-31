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

//- mobility computation 
surfaceTensorField Maf ("Maf",Kf*kraf/mua);
surfaceTensorField Laf ("Laf",rhoa*Kf*kraf/mua);	
surfaceTensorField Mbf ("Mbf",Kf*krbf/mub);
surfaceTensorField Lbf ("Lbf",rhob*Kf*krbf/mub);
surfaceTensorField Mf ("Mf",Maf+Mbf);
surfaceTensorField Lf ("Lf",Laf+Lbf);
surfaceScalarField Fbf ("Fbf",(krbf/mub) / ((kraf/mua) + (krbf/mub)));
volScalarField Fb ("Fb",(krb/mub) / ( (kra/mua) + (krb/mub) ));

//- capillarity computation (pc)
pcModel().correct();
surfaceVectorField gradpc("gradpc",fvc::interpolate(pcModel().dpcdS()*fvc::grad(Sb),"pc"));

//- fluxes depending on saturation
surfaceScalarField phiG("phiG",(Lf & g) & mesh.Sf());
surfaceScalarField phiPc("phiPc",(Mbf&gradpc) & mesh.Sf());
