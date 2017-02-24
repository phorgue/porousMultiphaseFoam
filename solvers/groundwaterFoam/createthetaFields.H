//- relative permeability (kr)
krModel->correct();
volScalarField krtheta = krModel->krb();
surfaceScalarField krthetaf ("krthetaf",fvc::interpolate(krtheta,"krtheta"));

//- mobility and fractional flow 
surfaceScalarField Mf ("Mf",rhotheta*mag(g)*Kf*krthetaf/mutheta);
surfaceScalarField Lf ("Lf",rhotheta*Kf*krthetaf/mutheta);

//- fluxes depending on saturation
surfaceScalarField phiG("phiG",(Lf * g) & mesh.Sf());

//- Test if gravity is present
if (mag(g).value() == 0)
{
    FatalErrorIn("createthetaFields.H")
        << " Magnitude of gravity mag(g) equal to zero " << abort(FatalError);
}

//- null Pc field for darcy velocity boundary conditions
surfaceScalarField phiPcNull("phiPc",0*phiG);
