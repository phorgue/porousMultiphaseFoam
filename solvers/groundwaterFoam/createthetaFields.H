//- relative permeability (kr)
krModel->correctb();
const volScalarField& krtheta = krModel->krb();
surfaceScalarField krthetaf ("krthetaf",fvc::interpolate(krtheta,"krtheta"));

//- mobility and fractional flow
surfaceScalarField Lf ("Lf",rhotheta*Kf*krthetaf/mutheta);
surfaceScalarField Mf ("Mf",mag(g)*Lf);

//- fluxes depending on saturation
surfaceScalarField phiG("phiG",(Lf * g) & mesh.Sf());
Utheta.correctBoundaryConditions();
forAll(mesh.boundary(),patchi)
{
    if (isA< fixedValueFvPatchField<vector> >(Utheta.boundaryField()[patchi]))
    {
        phi.boundaryFieldRef()[patchi] = Utheta.boundaryField()[patchi] & mesh.Sf().boundaryField()[patchi];
    }
}

//- Test if gravity is present
if (mag(g).value() == 0)
{
    FatalErrorIn("createthetaFields.H")
        << " Magnitude of gravity mag(g) equal to zero " << abort(FatalError);
}

//- null Pc field for darcy velocity boundary conditions
surfaceScalarField phiPcNull("phiPc",0*phiG);

//- h equation residual
volScalarField ResiduN = Ss*pcModel->Se() * fvc::ddt(h)
           + massConservativeTerms * fvc::ddt(theta)
           + (1 - massConservativeTerms) * pcModel->Ch() * fvc::ddt(h)
           - fvc::laplacian(Mf,h)
           + fvc::div(phiG)
           + sourceTerm
           ;
