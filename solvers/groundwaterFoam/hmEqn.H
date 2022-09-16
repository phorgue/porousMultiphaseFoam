{
    h_m.storePrevIter();

    fvScalarMatrix h_mEqn
        (
            //- transport terms
            - fvm::laplacian(Mf_m,h_m)
            + fvc::div(phiG_m)
            ==
            - sourceTerm_m
        );
    if (!steady)
    {
        h_mEqn +=
            //- accumulation terms
            (Ss_m*pcModel->Se() + pcModel_m->Ch()) * fvm::ddt(h_m)
            //-mass conservative terms
            + massConservativeTerms * (
                pcModel->Ch()*(h.oldTime()-h.prevIter())
                + ( theta - theta.oldTime()))
            /runTime.deltaT();
    }

    #include "updateForcing.H"

    hEqn.solve();

    deltah = h-h.prevIter();
    forAll(fixedPotentialIDList,celli) deltah[fixedPotentialIDList[celli]] = 0;
    deltahIter = gMax(mag(deltah.internalField())());

}
