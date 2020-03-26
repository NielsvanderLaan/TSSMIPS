#include "benders.h"

bool Benders::add_cut(BendersCut &cut, Master::Solution sol, double tol)
{
  bool improper_cut = d_master.add_cut(cut, sol, tol);

  if (not improper_cut)   
  {
    vector<double> coef_y(d_n2);
    d_pslp.add_cglp_rows(cut.d_beta.data(), 1 + cut.d_tau, coef_y.data(), cut.d_alpha);
    d_agg.add_rows(cut);
  }

  return improper_cut;
}