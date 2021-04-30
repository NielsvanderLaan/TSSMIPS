#include "aggregator.h"

BendersCut Aggregator::zk_cut(Master::Solution sol, Master &master, bool lap_cuts, bool affine, size_t maxRounds, double tol, double rho_tol)
{
  double *x = sol.xVals.data();
  double theta = sol.thetaVal;
  double rho = theta;
  double cRho = tol + 1;    // this choice ensures that the while loop is always entered

  BendersCut cut;
  size_t iter = 0;
  cerr << "fp iteration\n";
  size_t S = d_zk.size();
  while (cRho > rho_tol)
  {
    auto t1 = chrono::high_resolution_clock::now();
    cRho = -rho;
    cut = BendersCut{ 0, vector<double>(d_n1, 0.0), 0 };

    double lptime = 0;
    double avg_lptime = 0;    // per round
    double cuttime = 0;
    double avg_cuttime = 0;   // per round
    double nrounds = 0;
    double cpu_time = 0;

#pragma omp parallel for reduction(sum : cut) reduction(+:cRho, lptime, avg_lptime, cuttime, avg_cuttime, nrounds, cpu_time)
    for (size_t s = 0; s < d_zk.size(); ++s)
    {
      auto t1_loop = chrono::high_resolution_clock::now();
      double lptime_lc = 0;
      double cuttime_lc = 0;
      double nrounds_lc = 0;
      d_zk[s].update(x,affine ? GRB_INFINITY : rho);
      if (not d_zk[s].solve(x, theta, affine ? GRB_INFINITY : rho, master, maxRounds,
                            nrounds_lc, lptime_lc, cuttime_lc,
                            not lap_cuts))
      {
        exit(1561);
      }

      nrounds += nrounds_lc;
      lptime += lptime_lc;
      avg_lptime += lptime_lc / nrounds_lc;
      cuttime += cuttime_lc;
      avg_cuttime += cuttime_lc / nrounds_lc;

      cut += d_zk[s].subgradient() * d_probs[s];
      cRho += d_probs[s] * d_zk[s].d_objVal;

      auto t2_loop = chrono::high_resolution_clock::now();
      cpu_time += chrono::duration_cast<chrono::milliseconds>(t2_loop - t1_loop).count() / 1000.0;
    }

    if (affine)
      break;

    rho += cRho / (1 + cut.d_tau);

    auto t2 = chrono::high_resolution_clock::now();
    double time = chrono::duration_cast<chrono::milliseconds>(t2 - t1).count() / 1000.0;
    cerr << iter << " T: " << time <<
         " cpu_time: " << cpu_time <<
         " #cuts: "  << nrounds / S <<
         " T(LP): " << lptime / S <<
         " Ta(LP): " << avg_lptime / S <<
         " T(CGLP): "<< cuttime / S <<
         " Ta(CGLP): " << avg_cuttime / S << endl;
    ++iter;
  }

#pragma omp parallel for
  for (size_t s = 0; s < d_zk.size(); ++s)
    d_zk[s].clear();      // reduces memory usage

  return cut;
}

