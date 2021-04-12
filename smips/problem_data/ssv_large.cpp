#include "problem.h"

void Problem::ssv_large(size_t n1, size_t n2, size_t S, bool fs_continuous, bool ss_binary)
{
  d_n1 = n1; d_m1 = 0; d_fs_leq = 0; d_fs_geq = 0;
  d_p1 = fs_continuous? 0 : n1;
  d_n2 = n2; d_p2 = n2; d_m2 = n1; d_ss_leq = d_m2; d_ss_geq = 0;
  d_S = S;
  d_fix_rec = false;
  d_fix_tech = false;

  vector<double> l1(d_n1, 0.0);
  vector<double> l2(d_n2, 0.0);
  vector<double> u1(d_n1, 5.0);
  double ub = ss_binary ? 1.0 : GRB_INFINITY;
  vector<double> u2(d_n2, ub);

  vector<double> fs_costs(d_n1, 0.0);
  double step = 4.0 / (d_n1 - 1);
  for (size_t var = 0; var != d_n1; ++var)
    fs_costs[var] = -1 - step * var;


  vector<vector<double>> identity;
  for (size_t con = 0; con != d_m2; ++con)
  {
    vector<double> row(d_n1, 0.0);
    row[con] = 1.0;
    identity.push_back(row);
  }
  vector<vector<double>> tech;
  double frac = 1.0 / (d_n1 + 1);
  for (size_t con = 0; con != d_m2; ++con)
  {
    vector<double> row(d_n1, frac);
    row[con] =  2 * frac;
    tech.push_back(row);
  }

  for (size_t s = 0; s != d_S; ++s)
  {
    d_q_omega.push_back(d_gen.rand_unif_vec(d_n2, -30, -16));
    d_W_omega.push_back(d_gen.rand_unif_mat(d_m2, d_n2, 1, 6));
    d_omega.push_back(d_gen.unif_real_vec(d_m2, 5, 15));

    d_T_omega.push_back(d_gen.uni() > 0.5 ?
                        identity:
                        tech);
  }



  d_L = -450;

  d_l1 = l1;
  d_l2 = l2;
  d_u1 = u1;
  d_u2 = u2;

  d_c = fs_costs;
  d_q = d_q_omega.front();
  d_Wmat = d_W_omega.front();
  d_Tmat = d_T_omega.front();
  d_probs = vector<double> (d_S, 1.0 / d_S);

}