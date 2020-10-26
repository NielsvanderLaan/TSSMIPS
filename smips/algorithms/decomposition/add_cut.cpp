#include "benders.h"
#include <fstream>

bool Benders::add_cut(BendersCut &cut, Master::Solution sol, double tol, bool force)
{
  if (cut.d_tau > 0)
    cut.scale();

  bool improper_cut = not force;
  if (force)
    d_master.add_cut(cut);
  else
    improper_cut = d_master.add_cut(cut, sol, tol);

  if (not improper_cut)   
  {
    vector<double> coef_y(d_n2, 0.0);
    d_pslp.add_cglp_rows(cut.d_beta.data(), 1 + cut.d_tau, coef_y.data(), cut.d_alpha);
    d_agg.add_rows(cut);
    d_lr.add_cut(cut);
    d_def.add(cut);
  }

  return improper_cut;
}