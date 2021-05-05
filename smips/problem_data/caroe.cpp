#include "problem.h"

void Problem::caroe(size_t S)
{
  d_n1 = 1;
  d_m1 = 0;
  d_p1 = 0;
  d_fs_leq = 0;
  d_fs_geq = 0;

  d_n2 = 1;
  d_p2 = 1;
  d_m2 = 1;
  d_ss_leq = 0;
  d_ss_geq = 1;

  d_L = -2.0;

  d_S = S;
  d_probs = vector<double> (S, 1.0 / S);
  double step = (1.0 / 32.0) / (1.0 + S / 2.0);

  for (size_t s = 0; s != S / 2; ++s)
  {
    d_omega.push_back(vector<double> {0.25 - step * (s + 1)});
    d_omega.push_back(vector<double> {step * (s + 1)});
  }
  double low = 0;

  d_c.push_back(3.0);
  d_q.push_back(-2.0);

  d_Tmat.push_back(vector<double> { 1.0 });
  d_Wmat.push_back(vector<double> { -0.5 });
  d_l1.push_back(low);
  d_u1.push_back(1.0);
  d_l2.push_back(0.0);
  d_u2.push_back(1.0);

  /*
  d_L = 0;
  d_n2 = 2;
  d_q.push_back(2.0);
  d_l2.push_back(1.0);
  d_u2.push_back(1.0);
  d_Wmat[0].push_back(0.0);
  */


}