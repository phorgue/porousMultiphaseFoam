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
Fbf = (krbf/mub) / ( (kraf/mua) + (krbf/mub) );

//- for source term in saturation equation
Fb = (krb/mub) / ( (kra/mua) + (krb/mub) );

//- capillarity computation
pcModel().correct();
gradpc = fvc::interpolate(pcModel().dpcdS()*fvc::grad(Sb),"pc");

//- compute fluxes 
phiPc = (Mbf & gradpc) & mesh.Sf();
phiG = (Lf & g) & mesh.Sf();
