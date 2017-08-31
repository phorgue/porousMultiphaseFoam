//- relative permeability computation
krModel->correct();
kra = krModel->kra();
krb = krModel->krb();
dkradS = krModel->dkradS();
dkrbdS = krModel->dkrbdS();

kraf = fvc::interpolate(kra,"kra");
krbf = fvc::interpolate(krb,"krb");
dkrafdS = fvc::interpolate(dkradS,"kra");
dkrbfdS = fvc::interpolate(dkrbdS,"krb");

//- mobility computation 
Maf = Kf*kraf/mua;
Laf = rhoa*Kf*kraf/mua;	
Mbf = Kf*krbf/mub;
Lbf = rhob*Kf*krbf/mub;
Mf = Maf+Mbf;
Lf = Laf+Lbf;
Fbf = Mbf/Mf;

//- for source term in saturation equation
Fb = (krb/mub) / ( (kra/mua) + (krb/mub) );

//- compute fluxes
pcModel().correct();
phiPc = Mbf * fvc::interpolate(pcModel().dpcdS(),"pc")* fvc::snGrad(Sb) * mesh.magSf();
phiG = (Lf * g) & mesh.Sf();
