const bool subCycling = runTime.controlDict().lookupOrDefault("subCycling",false);

scalar dSmax(runTime.controlDict().lookupOrDefault<scalar>("dSmax",0.));
