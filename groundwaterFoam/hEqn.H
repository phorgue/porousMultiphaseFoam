{
 
    h.storePrevIter();

    fvScalarMatrix hEqn
        (
           pcModel->Ch() * fvm::ddt(h)
           - fvm::laplacian(Mf,h)
           + fvc::div(phiG)
        );

    hEqn.solve();
    
    resPicard = gMax((mag(h-h.prevIter()))->internalField());
    Info << "Picard iteration : " << iterPicard+1 << " residual = " << resPicard << endl; 

    phi = phiG-(Mf*fvc::snGrad(h))*mesh.magSf();
    Utheta = fvc::reconstruct(phi);


}
