//- update saturation 
theta = pcModel->correctAndSb(h);

//- relative permeability computation
krModel->correctb();
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

//- Compute h equation residual
ResiduN = Ss*pcModel->Se() * fvc::ddt(h)
    + massConservativeTerms * fvc::ddt(theta)
    + (1 - massConservativeTerms) * pcModel->Ch() * fvc::ddt(h)
    - fvc::laplacian(Mf,h)
    + fvc::div(phiG)
    + sourceTerm;
hEqnResidual = gMax(mag(ResiduN.internalField())());

