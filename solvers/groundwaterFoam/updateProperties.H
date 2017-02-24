//- update saturation 
theta == pcModel->correctAndSb(h);

//- relative permeability computation
krModel->correct();
krtheta = krModel->krb();
krthetaf = fvc::interpolate(krtheta,"krtheta");

//- mobility computation 
Mf = rhotheta*mag(g)*Kf*krthetaf/mutheta;
Lf = rhotheta*Kf*krthetaf/mutheta;

//- compute fluxes
phiG = (Lf * g) & mesh.Sf();
