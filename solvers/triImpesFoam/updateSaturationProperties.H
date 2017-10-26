//- relative permeability computation
krModel->correct();
kra = krModel->kra();
kro = krModel->kro();
krw = krModel->krw();

kraf = fvc::interpolate(kra,"kra");
krof = fvc::interpolate(kro,"kro");
krwf = fvc::interpolate(krw,"krw");

//- mobility computation 
Maf = Kf*kraf/mua;
Laf = rhoa*Kf*kraf/mua;
Mof = Kf*krof/muo;
Lof = rhoo*Kf*krof/muo;
Mwf = Kf*krwf/muw;
Lwf = rhow*Kf*krwf/muw;
Mf = Maf+Mwf+Mof;
Lf = Laf+Lwf+Lof;

//- compute fluxes
pcModel().correct();
phiPc_o = Mof * fvc::interpolate(pcModel().dpcaodS(),"pc")* fvc::snGrad(So) * mesh.magSf();
phiPc_w = Mwf * fvc::interpolate(pcModel().dpcaodS()+pcModel().dpcowdS(),"pc")* fvc::snGrad(Sw) * mesh.magSf();
//phiPc_w = Mwf * fvc::interpolate(pcModel().dpcawdS(),"pc")* fvc::snGrad(Sw) * mesh.magSf();
phiPc = phiPc_o + phiPc_w;
phiG = (Lf * g) & mesh.Sf();
