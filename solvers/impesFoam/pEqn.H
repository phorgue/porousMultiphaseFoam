{
    fvScalarMatrix pEqn
        (
            fvm::laplacian(-Mf, p) + fvc::div(phiG)
            // capillary term
            + fvc::div(phiPc)*activateCapillarity
            ==
            // event source terms
            - sourceTerm
        );

    pEqn.solve();

    phiP = pEqn.flux();

    phi = phiP+phiG+phiPc*activateCapillarity;

    phib == Fbf*phiP + (Lbf/Lf)*phiG + phiPc*activateCapillarity;
    phia == phi - phib;

    U = fvc::reconstruct(phi);
    U.correctBoundaryConditions();

    Ub = fvc::reconstruct(phib);
    Ua = U-Ub;

    Ub.correctBoundaryConditions();  
    Ua.correctBoundaryConditions();

}
