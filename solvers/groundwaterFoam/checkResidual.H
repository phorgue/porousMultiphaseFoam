{
    scalar hEqnResidualPrev = hEqnResidual;

    #include "updateProperties.H"

    if (hEqnResidual > 1)
    {
        iterPicard = maxIterPicard;
        iterNewton = maxIterNewton;
        Info << "Residual too large => reducting time step" << endl;
    }
    else if ( hEqnResidualPrev < hEqnResidual && iterNewton > 1)
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
            Info << "Non convergence of the algorithm to reduce large  residual = " << hEqnResidual << " => reducing time step" << endl;
        }
    }
}
