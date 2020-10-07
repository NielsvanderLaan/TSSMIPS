#include "fenchel.h"

Fenchel::Fenchel(const Fenchel &other)
:
  d_problem(other.d_problem),
  d_mp(other.d_mp),
  d_beta(other.d_beta.size()),
  d_sub(other.d_sub),
  d_xvars(other.d_xvars.size()),
  d_points(other.d_points),
  d_rcut(other.d_rcut)
{
  GRBVar *mp_vars = d_mp.getVars();
  d_alpha = mp_vars[0];
  copy_n(mp_vars + 1, d_beta.size(), d_beta.begin());
  d_kappa = mp_vars[d_beta.size() + 1];
  delete[] mp_vars;

  GRBVar *sub_vars = d_sub.getVars();
  copy_n(sub_vars, d_xvars.size(), d_xvars.begin());
  d_theta = sub_vars[d_beta.size()];
  delete[] sub_vars;

  d_mp.update();
  d_sub.update();
}