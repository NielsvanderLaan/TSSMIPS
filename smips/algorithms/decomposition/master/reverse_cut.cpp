#include "master.h"

void Master::reverse_cut(double UB)
{
  d_fenchel.reverse_cut(UB);
  if (d_rcut_idx == -1)
  {
    vector<double> minus_c = d_problem.d_c;
    for_each(minus_c.begin(), minus_c.end(), [](double &val){val *= -1;});
    BendersCut rcut {-UB, minus_c ,-2};
    add_cut(rcut);
    int nCons;
    GRBgetintattr(d_cmodel, "NumConstrs", &nCons);
    d_rcut_idx = nCons - 1;
  } else
    GRBsetdblattrelement(d_cmodel, "RHS", d_rcut_idx, -UB + d_L);

  GRBoptimize(d_cmodel);
}