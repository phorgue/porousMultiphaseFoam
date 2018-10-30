const dictionary& PicardDict = mesh.solutionDict().subDict("Picard");
const scalar tolPicard = PicardDict.lookupOrDefault<scalar>("tolerance", 1.);
scalar resPicard=GREAT;
scalar sousrelaxfactor(PicardDict.lookupOrDefault<scalar>("sousrelaxfactor",0.001));
