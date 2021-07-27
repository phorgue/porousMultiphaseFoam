{
    scalar hEqnResidualPrev = hEqnResidual;

    label iter_relax = 1;
    h.relax();
    #include "updateProperties.H"

    if ( hEqnResidualPrev < hEqnResidual )
    {
        iter_relax++;
        h.relax(0.5);
        Info << "Additional under-relaxation of the solution to reduce residual  " << endl;
        #include "updateProperties.H"
    }
    Info << "Head pressure equation residual = " << hEqnResidual << endl;
}
