{

    fvScalarMatrix pEqn
        (
            fvm::laplacian(-Mf, p) + fvc::div(phiG)
            // capillary term
            + fvc::div(phiPc)*activateCapillarity
            // wellbore terms
            - (-SrcExt*Wext+SrcInj*Winj)
        );

    pEqn.solve();

    phiP = pEqn.flux();

    phi = phiP+phiG+phiPc*activateCapillarity;

    U = fvc::reconstruct(phi);
    U.correctBoundaryConditions();

    phib == Fbf*phi+(((Lbf-Fbf*Lf)&g) & mesh.Sf())+(1-Fbf)*phiPc*activateCapillarity;
    phia == phi - phib;

    Ub = fvc::reconstruct(phib);
    Ua = U-Ub;

    Ua.correctBoundaryConditions();
    Ub.correctBoundaryConditions();  

}
