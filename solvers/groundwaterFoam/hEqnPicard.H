{
    porousModel->correct(h, steady, massConservative);

    h.storePrevIter();

    fvScalarMatrix hEqn
        (
            //- transport terms
            - fvm::laplacian(Mf,h)
            + fvc::div(phiG)
            ==
            - porousModel->source()
            - sourceTerm
        );
    if (!steady)
    {
        //- accumulation terms
        hEqn += (Ss*pcModel->Se() + pcModel->Ch()) * fvm::ddt(h);

        if (massConservative)
        {
            //-mass conservative terms
            hEqn += (pcModel->Ch()*(h.oldTime()-h.prevIter())
            + (theta - theta.oldTime())) / runTime.deltaT();
        }
    }

    #include "updateForcing.H"

    hEqn.solve();

    deltah = h-h.prevIter();
    forAll(fixedPotentialIDList,celli) deltah[fixedPotentialIDList[celli]] = 0;
    deltahIter = gMax(mag(deltah.internalField())());

}
