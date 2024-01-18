{
    porousModel->correct(h, steady, massConservative);

    h.storePrevIter();

    fvScalarMatrix hEqn
        (
            //- transport terms
            - fvm::laplacian(Mf,h)
            + fvc::div(phiG)
            ==
            - porousModel->exchangeTerm()
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

    if (fixedPotentialIDList.size() > 0) hEqn.setValues(fixedPotentialIDList,fixedPotentialValueList);

    hEqnResidual = hEqn.solve().initialResidual();
    if (hEqnResidual > tolerancePicard) h.relax();

    deltah = h-h.prevIter();
    forAll(fixedPotentialIDList,celli) deltah[fixedPotentialIDList[celli]] = 0;
    deltahIter = gMax(mag(deltah.internalField())());

}
