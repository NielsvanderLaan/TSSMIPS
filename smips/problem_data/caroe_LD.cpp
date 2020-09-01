#include "problem.h"

void Problem::caroe_LD(size_t S)
{
  d_fix_rec = false;

  d_n1 = 1;
  d_m1 = 0;
  d_p1 = 0;
  d_fs_leq = 0;
  d_fs_geq = 0;

  d_n2 = 1;
  d_p2 = 0;
  d_m2 = 2;
  d_ss_leq = 0;
  d_ss_geq = 2;

  d_L = -2;

  d_S = S;
  d_probs = vector<double> (S, 1.0 / S);
  double step = (1.0 / 32.0) / (1.0 + S / 2.0);

  double low = 0.25 - step;

  for (size_t s = 0; s != S / 2; ++s)
  {
    double h = 0.25 - step * (s + 1);
    d_omega.push_back(vector<double> {h, low});
    d_q_omega.push_back(vector<double> {-2.0});

    vector<vector<double>> Wmat;
    Wmat.push_back(vector<double> {-0.5});
    Wmat.push_back(vector<double> {low - 0.5 - h});
    d_W_omega.push_back(Wmat);


    h = step * (s + 1);
    d_omega.push_back(vector<double> {h, low});
    d_q_omega.push_back(vector<double> {-2.0});
    Wmat[1][0] = low - 0.5 - h;
    d_W_omega.push_back(Wmat);
  }

  d_c.push_back(3.0);

  d_Tmat.push_back(vector<double> { 1.0 });
  d_Tmat.push_back(vector<double> { 1.0 });
  d_Wmat.push_back(vector<double> { -0.5 });
  d_l1.push_back(low);
  d_u1.push_back(1.0);
  d_l2.push_back(0.0);
  d_u2.push_back(1.0);


}