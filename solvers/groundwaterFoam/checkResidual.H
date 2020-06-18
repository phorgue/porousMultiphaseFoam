{
    scalar hEqnResidualPrev = hEqnResidual;
    scalar hEqnResidualSignedPrev = hEqnResidualSigned;

    #include "updateProperties.H"

    if ( hEqnResidualPrev < hEqnResidual && (hEqnResidualSignedPrev * hEqnResidualSigned) > 0 )
    {
        label relaxation_iterator = 0;
        while (relaxation_iterator < 3 && hEqnResidual > hEqnResidualPrev)
        {
            relaxation_iterator++;
            Info << "Under-relaxation of the solution to reduce large residual = " << hEqnResidual << endl;
            h.relax(0.5);
        #include "updateProperties.H"
        }
        if (relaxation_iterator == 3)
        {
            iterNewton = maxIterNewton;
            iterPicard = maxIterPicard;
            Info << "Non convergence of the algorithm => reducing time step" << endl;
        }
    }

}
