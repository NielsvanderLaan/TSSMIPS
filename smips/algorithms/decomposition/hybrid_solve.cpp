#include "benders.h"

Benders::Bounds Benders::hybrid_solve(vector<Type> types, bool force_int, size_t max_rounds,
                                      double upper_bound, double tol, double time_limit)
{
  size_t gmi_cuts = 0;
  vector<size_t> nCuts(types.size(), 0);
  size_t round = 0;       // rounds of gmi cuts added

  double LB;
  bool branch = false;

  auto t1 = chrono::high_resolution_clock::now();
  while (true)
  {
    start:
    auto t2 = chrono::high_resolution_clock::now();
    if (chrono::duration_cast<chrono::milliseconds>(t2 - t1).count() / 1000.0 > time_limit)
    {
      cout << "OOT\n";
      break;
    }

    Master::Solution sol = d_master.solve(tol);
    if (sol.infeasible)
    {
      cout << "mp infeasible" << endl;
      return Bounds{GRB_INFINITY, GRB_INFINITY};
    }

    vector<double> x = sol.xVals;

    LB = get_lb();

    if (LB > upper_bound)// || LB > UB - 1e-8)
    {
      cout << "LB > upper_bound (LB = " << LB << ")" << endl;
      break;
    }
    bool int_feas = all_of(x.begin(), x.begin() + d_p1, [](double val){ return is_integer(val); });

    if (not int_feas && round < max_rounds)
    {
      size_t nCuts = round_of_cuts(sol, 1e-4);
      gmi_cuts += nCuts;
      if (nCuts > 0)      // at least one cut was added
      {
        print("added " << nCuts <<  " gmi cuts\n");
        ++round;
        continue;
      }
    }

    if (not int_feas)
    {
      copy(x.begin(), x.end(), d_xvals);
      branch = true;
      if (force_int)  break;
    }

    double cx = inner_product(d_problem.d_c.data(), d_problem.d_c.data() + d_n1, x.begin(), 0.0);
    vector<double> vx = d_agg.compute_vwx(x.data());
    double Qx = accumulate(vx.begin(), vx.end(),0.0) / d_S;
    if (int_feas && cx + Qx < d_UB)
    {
      copy(x.begin(), x.end(), d_incumbent);
      update(cx + Qx);    // updates UB and calls reverse_cut
    }
    print("LB: " << LB << ". UB: " << d_UB << endl);

    for (size_t idx = 0; idx != types.size(); ++idx)
    {
      BendersCut cut = compute_cut(types[idx], sol, int_feas, vx, tol);
      if (not add_cut(cut, sol, tol))
      {
        print("added " << name(types[idx]) << '\n');
        ++nCuts[idx];
        goto start;
      }
    }

    copy(x.begin(), x.end(), d_xvals);
    cout << "no improvement possible\n";
    branch = true;
    break;
  }
  cout << "gmi cuts: " << gmi_cuts << " (" << round << " rounds)\n";
  for (size_t idx = 0; idx != types.size(); ++idx)
    cout << name(types[idx]) << "s: " << nCuts[idx] << '\n';

  cout << "LB: " << LB << ". UB: " << d_UB << endl;
  return Bounds { LB, d_UB, branch };
}





