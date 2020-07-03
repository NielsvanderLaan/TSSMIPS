#include "master.h"

vector<double> Master::extract_row(size_t row)
{
  int nVars = 1 + d_n1 + d_nSlacks;
  
  int nCons;
  GRBgetintattr(d_cmodel, "NumConstrs", &nCons);

  int len = nVars + nCons;
  int inds[len];
  double vals[len];
  GRBsvec tab_row{ len, inds, vals };
  
  GRBBinvRowi(d_cmodel, row, &tab_row);

        // converting to std::vector
  vector<double> tab_row_vec(nVars);
  for (size_t nz = 0; nz != tab_row.len; ++nz)
  {
    size_t idx = tab_row.ind[nz];
    if (idx < nVars)
      tab_row_vec[idx] = tab_row.val[nz];
  }
  
  return tab_row_vec;
}