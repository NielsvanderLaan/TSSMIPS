#include "lagrangian.h"

BendersCut Lagrangian::lr_cut(size_t s, double *x, double vwx, vector<double> &pi)
{
  update(s, pi);
  double Lpiw = solve();

  double step;
  double m = 1.0;
  //cout << "start:\n";
  //for (size_t k = 1; k != 101; ++k)
  size_t k = 0;
  while (true)
  {
    ++k;
    //cout << "k = " << k << ": ";
    //for_each(pi.begin(), pi.end(), [](double val){cout << val << ' ';});
    //cout << endl;
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
    step = (1.0 + m) / (k + m) * (vwx - qmu_k) / gk_squared;
    if (step * step * gk_squared < 1e-4)
      break;

    //cout << "step = " << step << '\n';

    for (size_t var = 0; var != d_n1; ++var)
      pi[var] += step * gk[var];
    update_pi(pi);
    Lpiw = solve();
  }

  return BendersCut {Lpiw, pi, 0.0};
}