#include "benders.h"

Benders::Bounds Benders::hybrid_solve(vector<Type> types, bool force_int, int max_rounds,
                                      double upper_bound, double tol, double time_limit,
                                      bool rcuts, bool fenchel)
{
  size_t gmi_cuts = 0;
  size_t fenchel_cuts = 0;
  double gmi_time = 0.0;
  double fenchel_time = 0.0;
  vector<size_t> nCuts(types.size(), 0);
  vector<double> times(types.size(), 0.0);
  size_t evaluations = 0;
  double eval_time = 0.0;
  size_t round = 0;       // rounds of gmi cuts added

  double LB = -GRB_INFINITY;
  bool branch = false;

  size_t nStall = 10;
  double stall_tol = 1e-4;      // relative (if LB does not improve by 1e-4*100% for nStall iterations, then skip hierarchy)
  list<double> recent_lbs;

  auto t1 = chrono::high_resolution_clock::now();
  while (true)
  {
    start:
    if (chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - t1).count() / 1000.0 > time_limit)
    {
      cout << "OOT\n";
      break;
    }

    Master::Solution sol = d_master.solve(tol);
    if (sol.infeasible)
    {
      cout << "mp infeasible" << endl;
      LB = GRB_INFINITY;
      branch = false;
      break;
    }

    vector<double> x = sol.xVals;
    LB = get_lb();

    if (LB > upper_bound)
    {
      cout << "LB > upper_bound" << endl;
      break;
    }

    bool int_feas = all_of(x.begin(), x.begin() + d_p1, [](double val){ return is_integer(val); });


    if (not int_feas)
    {
      if (fenchel && fenchel_cuts != max_rounds)
      {
        auto before = chrono::high_resolution_clock::now();
        BendersCut cut = d_master.fenchel_cut(sol, tol);
        double comp_time = chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - before).count() / 1000.0;
        print("computed Fenchel cut (" << comp_time << "s)\n");
        if (not add_cut(cut, sol, tol))
        {
          print("added Fenchel cut\n");
          fenchel_time += comp_time;
          ++fenchel_cuts;
          continue;
        }
      }
      if (not fenchel and round != max_rounds)
      {
        auto before = chrono::high_resolution_clock::now();
        size_t nCuts = round_of_cuts(sol, 1e-4);
        double time = chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - before).count() / 1000.0;
        print("computed round of GMI cuts (" << time << "s)\n");
        if (nCuts > 0)      // at least one cut was added
        {
          print("added " << nCuts <<  " gmi cuts\n");
          gmi_time += time;
          gmi_cuts += nCuts;
          ++round;
          continue;
        }
      }
    }

    if (not int_feas and force_int)
    {
      copy(x.begin(), x.end(), d_xvals);
      branch = true;
      cout << "not integer feasible. #rounds of cuts: " << round << ". #fenchel cuts: " << fenchel_cuts << '\n';
      break;
    }

    auto before = chrono::high_resolution_clock::now();
    double cx = inner_product(d_problem.d_c.data(), d_problem.d_c.data() + d_n1, x.begin(), 0.0);
    vector<double> vx = d_agg.compute_vwx(x.data());
    double Qx = accumulate(vx.begin(), vx.end(),0.0) / d_S;
    double time = chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - before).count() / 1000.0;
    ++evaluations;
    print("evaluated cx + Q(x) (" << time << "s)\n");
    eval_time += time;

    if (int_feas && cx + Qx < d_UB)
    {
      copy(x.begin(), x.end(), d_incumbent);
      update(cx + Qx, rcuts);    // updates UB and calls reverse_cut
    }
    print("LB: " << LB << ". UB: " << d_UB << endl);

    recent_lbs.push_front(LB);
    if (recent_lbs.size() > nStall + 1) recent_lbs.pop_back();
    double old_LB = recent_lbs.back();
    size_t start = recent_lbs.size() > nStall && LB - old_LB < stall_tol * abs(old_LB)  ? types.size() - 1 : 0;
    if (start != 0) cout << "BREAKING PRIORITY RULES\n";

    for (size_t idx = start; idx != types.size(); ++idx)
    {
      auto before = chrono::high_resolution_clock::now();
      BendersCut cut = compute_cut(types[idx], sol, int_feas, vx, tol);
      double time = chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - before).count() / 1000.0;
      print("computed " << name(types[idx]) << " (" << time << "s)" << '\n');

      cout << "tau = " << cut.d_tau << '\n';
      if (not add_cut(cut, sol, tol))
      {
        print("added " << name(types[idx]) << '\n');
        ++nCuts[idx];
        times[idx] += time;
        goto start;
      }
    }

    copy(x.begin(), x.end(), d_xvals);
    cout << "no improvement possible\n";
    branch = true;
    break;
  }

  auto t2 = chrono::high_resolution_clock::now();
  cout << "LB: " << LB << ". UB: " << d_UB << '\n';
  cout << "computation time: " << chrono::duration_cast<chrono::milliseconds>(t2 - t1).count() / 1000.0 << '\n';
  cout << "evaluations: " << evaluations << " (" << avg(eval_time, evaluations) << "s)\n";
  cout << "GMI cuts: " << gmi_cuts << " (" << round << " rounds, " << avg(gmi_time, round) << "s)\n";
  cout << "Fenchel cuts: " << fenchel_cuts << " (" << avg(fenchel_time, fenchel_cuts) << "s)\n";
  for (size_t idx = 0; idx != types.size(); ++idx)
    cout << name(types[idx]) << "s: " << nCuts[idx] << " (" << avg(times[idx], nCuts[idx]) << "s)\n";

  return Bounds { LB, d_UB, branch };
}







