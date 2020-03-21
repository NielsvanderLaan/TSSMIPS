#include "master.h"

vector<BendersCut> Master::round_of_cuts()
{
  int nCons;
  GRBgetintattr(d_cmodel, "NumConstrs", &nCons);
  int bhead[nCons];
  GRBgetBasisHead(d_cmodel, bhead);        // extract basis info  
  
  double x[d_p1];
  GRBgetdblattrarray(d_cmodel, "X", 0, d_p1 + 1, x);   
  
  vector<BendersCut> cuts;
  
  for (size_t row = 0; row != nCons; ++row)    // loop over rows of simplex tableau 
  {
    int var = bhead[row];                // index of corresponding basic variable
    if (var == 0 || var > d_p1)                     // check if variable has to be integer
      continue;
      
    double val = x[var]; 
    if (is_integer(val))                   // check if integer requirement is violated
      continue;                            // then do not derive a cut
        
    cuts.push_back(gmi_cut(row, val));
  }

  return cuts;
}