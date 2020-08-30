#include "lagrangian.h"

BendersCut Lagrangian::lr_cut(size_t s, double *x, double vwx, vector<double> &pi)
{
  update(s, pi);
  double Lpiw = solve();

  double step;
  int m = 1;
  double previous = 0.0;

  for (size_t k = 1; k != 101; ++k)
  {
    vector<double> gk = z_vals();     // subgradient
    double gk_squared = 0.0;          // squared L2 norm of gk
    for (size_t var = 0; var != gk.size(); ++var)
    {
      gk[var] -= x[var];
      gk_squared += gk[var] * gk[var];
    }
    if (gk_squared < 1e-4)
      break;
    double qmu_k = -inner_product(pi.begin(), pi.end(), x, -Lpiw);
    if (abs(qmu_k - previous) < 1e-4 and k > 1)
      break;
    previous = qmu_k;

    step = (1.0 + m) / (k + m) * (vwx - qmu_k) / gk_squared;

    for (size_t var = 0; var != d_n1; ++var)
      pi[var] += step * gk[var];
    update_pi(pi);
    Lpiw = solve();
  }

  return BendersCut {Lpiw, pi, 0.0};
}