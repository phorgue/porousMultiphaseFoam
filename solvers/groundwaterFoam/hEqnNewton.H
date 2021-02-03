{
    h.storePrevIter();
    const volScalarField& dkrthetadS = krModel ->dkrbdS();
    volScalarField dLdS = pcModel->Ch() * rhotheta * K * dkrthetadS / mutheta;
    volScalarField dMdS = mag(g) * dLdS;

    fvScalarMatrix deltahEqn_hGrad(dMdS  * fvm::div(fvc::snGrad(h) * mesh.magSf(),deltah));
    fvScalarMatrix deltahEqn_grav( dLdS * fvm::div(-g & mesh.Sf(),deltah));

    deltahEqn_hGrad.diag() *= -1;
    deltahEqn_grav.diag() *= -1;

    fvScalarMatrix deltahEqn
        (
            fvm::Sp(
                (
                    (Ss*pcModel->Ch()) * (h - h.oldTime())
                    + Ss*pcModel->Se()
                    + pcModel->Ch()
                )/runTime.deltaT()
                ,deltah)
            - fvm::laplacian(Mf,deltah)
            + deltahEqn_hGrad
            + deltahEqn_grav
            ==
            - ResiduN
        );

    scalarField forInversion = deltahEqn.upper();
    deltahEqn.upper() = deltahEqn.lower();
    deltahEqn.lower() = forInversion;

    forAll(mesh.boundary(),patchi)
    {
        if (deltah.boundaryField().types()[patchi] == "darcyGradPressure")
        {
            deltahEqn.internalCoeffs()[patchi] = 0;
            deltahEqn.boundaryCoeffs()[patchi] = 0;
        }
    }

    deltahEqn.solve();

    deltahIter = max(mag(deltah.internalField())).value();
    h = h.prevIter() + deltah;
    h.correctBoundaryConditions();

}
