//- update saturation 
theta = pcModel->correctAndSb(h);

//- relative permeability computation
krModel->correctkrb(theta);
krthetaf = fvc::interpolate(krtheta,"krtheta");

//- mobility computation 
Lf = rhotheta*Kf*krthetaf/mutheta;
Mf = mag(g)*Lf;

//- compute fluxes
phiG = (Lf * g) & mesh.Sf();
phi = phiG-(Mf*fvc::snGrad(h))*mesh.magSf();
Utheta = fvc::reconstruct(phi);
Utheta.correctBoundaryConditions();
forAll(mesh.boundary(),patchi)
{
    if (isA< fixedValueFvPatchField<vector> >(Utheta.boundaryField()[patchi]))
    {
        phi.boundaryFieldRef()[patchi] = Utheta.boundaryField()[patchi] & mesh.Sf().boundaryField()[patchi];
    }
}
h.correctBoundaryConditions();

//- Compute h equation residual
ResiduN = - fvc::laplacian(Mf,h) + fvc::div(phiG) + porousModel->source() + sourceTerm;
if (!steady)
{
    ResiduN += Ss*pcModel->Se() * fvc::ddt(h);
    if (massConservative) ResiduN += fvc::ddt(theta);
    else ResiduN += pcModel->Ch() * fvc::ddt(h);
}

forAll(fixedPotentialIDList, pointi) ResiduN[fixedPotentialIDList[pointi]]=0;
hEqnResidual = gMax(mag(ResiduN.internalField())());
