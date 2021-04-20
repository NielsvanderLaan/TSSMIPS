#include "problem.h"

void Problem::ssv_large(int n1, int n2, int m2, size_t S, bool fs_continuous, bool ss_binary)
{
  d_n1 = n1; d_m1 = 0; d_fs_leq = 0; d_fs_geq = 0;
  d_p1 = fs_continuous? 0 : n1;
  d_n2 = n2; d_p2 = n2; d_m2 = m2; d_ss_leq = d_m2; d_ss_geq = 0;
  d_S = S;
  d_fix_rec = false;
  d_fix_tech = false;

  vector<double> l1(d_n1, 0.0);
  vector<double> l2(d_n2, 0.0);
  vector<double> u1(d_n1, 5.0);
  double ub = ss_binary ? 1.0 : GRB_INFINITY;
  vector<double> u2(d_n2, ub);

  vector<double> fs_costs(d_n1);
  double step = 4.0 / (d_n1 - 1);
  for (size_t var = 0; var != d_n1; ++var)
    fs_costs[var] = -1 - step * var;


  vector<vector<double>> identity, tech;
  if (n1 > m2)
  {
    int nz = n1 -m2+ 1;
    double frac = 1.0 / (n1 + nz);
    for (size_t con = 0; con != m2; ++con)
    {
      vector<double> row(n1, 0.0);
      fill_n(row.begin() + con, nz, 1.0 / nz);
      identity.push_back(row);
    }

    for (size_t con = 0; con != m2; ++con)
    {
      vector<double> row(n1, frac);
      fill_n(row.begin() + con, nz, 2 * frac);
      tech.push_back(row);
    }
  } else
  {
    int nz = min(m2 - n1 + 1, n1);
    int ndiags = 0;
    for (int con = 0; con != nz; ++con)
    {
      ++ndiags;
      vector<double> row(n1, 0.0);
      fill_n(row.begin(), ndiags, 1.0 / ndiags);
      identity.push_back(row);

      vector<double> t_row(n1, 1.0 / (n1 + ndiags));
      fill_n(t_row.begin(), ndiags, 2.0 / (n1 + ndiags));
      tech.push_back(t_row);
    }
    int nrows = m2 - 2 * nz;
    for (int con = 0; con < nrows; ++con)
    {
      vector<double> row(n1, 0.0);
      if (nz < n1)
        fill_n(row.begin() + con + 1, nz, 1.0 / nz);
      else
        fill_n(row.begin(), nz, 1.0 / nz);
      identity.push_back(row);

      vector<double> t_row(n1, 1.0/(n1 + nz));
      if (nz < n1)
        fill_n(t_row.begin() + con + 1, nz, 2.0 / (n1 + nz));
      else
        fill_n(t_row.begin(), nz, 2.0 / (n1 + nz));
      tech.push_back(t_row);

    }
    for (int con = nz; con != 0; --con)
    {
      vector<double> row(n1, 0.0);
      fill(row.end() - con, row.end(), 1.0 / ndiags);
      identity.push_back(row);

      vector<double> t_row(n1, 1.0 / (n1 + ndiags));
      fill(t_row.end() - con, t_row.end(), 2.0 / (n1 + ndiags));
      tech.push_back(t_row);
      --ndiags;
    }
  }

  int qlow = -(10*n1 + 10);
  int qhigh = -(10*n1 - 10);
  for (size_t s = 0; s != d_S; ++s)
  {
    d_q_omega.push_back(d_gen.rand_unif_vec(d_n2, qlow, qhigh));
    d_W_omega.push_back(d_gen.rand_unif_mat(d_m2, d_n2, 1, 6));
    d_omega.push_back(d_gen.unif_real_vec(d_m2, 5, 15));

    d_T_omega.push_back(d_gen.uni() > 0.5 ?
                        identity:
                        tech);
  }


  d_L = 15*qlow;

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