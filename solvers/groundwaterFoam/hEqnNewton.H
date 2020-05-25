{
    h.storePrevIter();
    
    const volScalarField& dkrthetadS = krModel ->dkrbdS();
    volScalarField LdS = pcModel->Ch() * rhotheta * K * dkrthetadS / mutheta;
    volScalarField MdS = mag(g) * LdS;

    fvScalarMatrix deltahEqn_hGrad(MdS  * fvm::div(fvc::snGrad(h) * mesh.magSf(),deltah));
    deltahEqn_hGrad.diag() *= -1;

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
            + LdS * fvm::div(-g & mesh.Sf(),deltah)
            ==
            - ResiduN
        );

    scalarField forInversion = deltahEqn.upper();
    deltahEqn.upper() = deltahEqn.lower();
    deltahEqn.lower() = forInversion;

    deltahEqn.solve();

    deltahIter = max(mag(deltah.internalField())).value();
    h = h.prevIter() + deltah;
    h.correctBoundaryConditions();

}
