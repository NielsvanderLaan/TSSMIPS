#include "aggregator.h"

BendersCut Aggregator::strong_cut(Master::Solution sol, vector<double> &vx, bool affine, double tol, bool int_feas, double rho_tol)
{
  double rho = sol.thetaVal;
  double *x = sol.xVals.data();
  double cRho = rho_tol + 1;
  BendersCut cut;

  double gap;
  bool first_time = true;
  cerr << "fp iteration\n";
  size_t iter = 0;
  int S = d_cgmips.size();
  while (cRho > rho_tol)
  {
    auto t1 = chrono::high_resolution_clock::now();

    cRho = -rho;
    cut = BendersCut{ 0, vector<double>(d_n1, 0.0), 0};
    gap = 0;

    double npoints = 0;
    double sptime = 0;
    double mptime = 0;
    double niter = 0;
    double avg_sptime = 0;
    double avg_mptime = 0;
#pragma omp parallel for reduction(sum : cut) reduction(+:cRho, gap, npoints, sptime, mptime, niter, avg_sptime, avg_mptime)
    for (size_t s = 0; s < S; ++s)
    {
      double prob = d_probs[s];
      double niter_lc, sptime_lc, mptime_lc;
      cut += d_cgmips[s].generate_cut(x, rho, first_time, vx[s], affine, tol, int_feas, gap,
                                      niter_lc, sptime_lc, mptime_lc) * prob;
      niter += niter_lc;
      mptime += mptime_lc;
      avg_mptime += mptime_lc / niter_lc;
      sptime += sptime_lc;
      avg_sptime += sptime_lc / niter_lc;
      npoints += d_cgmips[s].d_points.size();

      cRho -= prob * d_cgmips[s].mp_val();
    }


    gap /= d_cgmips.size();

    auto t2 = chrono::high_resolution_clock::now();
    double time = chrono::duration_cast<chrono::milliseconds>(t2 - t1).count() / 1000.0;
    cerr << iter << " T: " << time <<
         " #EP: " << npoints / S <<
         " #iter: "  << niter / S <<
         " T(CGMP): " << mptime / S <<
         " Ta(CGMP): " << avg_mptime / S <<
         " T(CGSP): "<< sptime / S <<
         " Ta(CGSP): " << avg_sptime / S << endl;

    if (affine)
      break;

    rho += cRho / (1 + cut.d_tau);
    first_time = false;
    ++iter;
  }

  for (size_t s = 0; s < d_cgmips.size(); ++s)
  {
    d_cgmips[s].update_mp();
    //d_cgmips[s].clear_mp();
  }


  if (gap > tol)
    cout << "strong_cut() gap: " << gap << '\n';

  return cut;
}
