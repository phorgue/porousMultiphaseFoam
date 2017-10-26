{
    fvScalarMatrix pEqn
        (
            fvm::laplacian(-Mf, p) + fvc::div(phiG)
            // capillary term
            + fvc::div(phiPc_o) * activateCapillarity
            + fvc::div(phiPc_w) * activateCapillarity          
            // wellbore terms
            - (SrcInjo*Winjo-Fo*SrcExto*Wexto)
            // - (-SrcExt*Wext+SrcInj*Winj)*activateWellbores 
        );

    pEqn.solve();

    phiP = pEqn.flux();

    phi = phiP+phiG+(phiPc_o+phiPc_w)*activateCapillarity;

    phio == (Mof/Mf)*phiP + (Lof/Lf)*phiG + phiPc_o*activateCapillarity;  
    phiw == (Mwf/Mf)*phiP + (Lwf/Lf)*phiG + phiPc_w*activateCapillarity;
    phia == phi - phiw - phio;

    U = fvc::reconstruct(phi);
    U.correctBoundaryConditions();

    Uw = fvc::reconstruct(phiw);
    Uo = fvc::reconstruct(phio);   
    Ua = U-Uw-Uo;

}
