{
    
    fvScalarMatrix hwaterEqn
        (
            fvm::ddt(hwater)
            + fvm::div(phi,hwater,"div(phi,hwater)") 
        );

    hwaterEqn.solve();
 
    dhdTmax = max(mag(fvc::ddt(hwater))).value();
    
    //- updating potential
    potential == hwater + z0;

    //- updating flux/velocity field
    phi == (-Mf * fvc::snGrad(potential)) * mesh.magSf();
    U = fvc::reconstruct(phi);

}
