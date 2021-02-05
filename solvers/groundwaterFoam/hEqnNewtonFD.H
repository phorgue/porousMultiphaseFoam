{
    h.storePrevIter();
    scalar epsilon = 1e-5;
    forAll(h,celli)
    {
        //- small variation of pressure head
        h[celli] -= epsilon;

        //- update moisture/Ch for considered cell
        scalar theta_backup = theta[celli];
        scalar Ch_backup = pcModel->Ch(celli);
        theta[celli] = pcModel->correctAndSb(h,celli);

        //- update relative permeability for considered cell
        scalar kr_backup = krModel->krb(celli);
        krModel->correctb(celli);

        krthetaf = fvc::interpolate(krtheta,"krtheta");
        Lf = rhotheta*Kf*krthetaf/mutheta;
        Mf = mag(g)*Lf;
        phiG = (Lf * g) & mesh.Sf();

        //- Compute new F_epsilon value
        volScalarField Feps = Ss*pcModel->Se() * fvc::ddt(h)
            + massConservativeTerms * fvc::ddt(theta)
            + (1 - massConservativeTerms) * pcModel->Ch() * fvc::ddt(h)
            - fvc::laplacian(Mf,h)
            + fvc::div(phiG)
            + sourceTerm;
        volScalarField dF = (ResiduN - Feps) / epsilon;

        //- store derivative values
        jacobian.storeColumn(dF,celli);

        //- restore values
        pcModel->setCh(celli,Ch_backup);
        krModel->setKrb(celli,kr_backup);
        theta[celli] = theta_backup;
        h[celli] += epsilon;

    }

    jacobian.matrix().source() = -ResiduN;
    jacobian.matrix().solve();

    deltahIter = max(mag(deltah.internalField())).value();
    h = h.prevIter() + deltah;
}
