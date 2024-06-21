//- Compute h equation residual
ResiduN = - fvc::laplacian(Mf,h) + fvc::div(phiG) + porousModel->exchangeTerm() + sourceTerm;
if (!steady)
{
    ResiduN += Ss*pcModel->Se() * fvc::ddt(h);
    if (massConservative) ResiduN += fvc::ddt(theta);
    else ResiduN += pcModel->Ch() * fvc::ddt(h);
}

forAll(fixedPotentialIDList, pointi) ResiduN[fixedPotentialIDList[pointi]]=0;
hEqnResidualMax = gMax(mag(ResiduN.internalField())());
