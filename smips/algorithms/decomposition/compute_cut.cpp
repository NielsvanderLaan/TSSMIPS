#include"benders.h"

BendersCut Benders::compute_cut(Type type, Master::Solution &sol, bool int_feas, vector<double> &vx, double tol, vector<double> alpha)
{
  switch (type)
  {
    case LP:        return d_agg.lp_cut(sol.xVals);
    case SB:        return d_agg.sb_cut(sol.xVals);
    case LR:        return d_agg.strong_cut(sol, vx, true, tol, int_feas);
    case LR_LAP:    return d_agg.zk_cut(sol, d_master, true, true);
    case SC_ZK:     return d_agg.zk_cut(sol, d_master, false, false);
    case SC_LAP:    return d_agg.zk_cut(sol, d_master, true, false);
    case SC_BAB:    return d_agg.bac_cut(sol, d_master, tol);
    case SC_RG:     return d_agg.strong_cut(sol, vx, false, tol, int_feas);
    case LBDA:
      if (alpha.size() != d_m2)
      {
        cerr << "alpha has wrong size\n";
        exit(1);
      }
      return lbdaCut(sol.xVals, alpha);
    default:
      cerr << "cut type unknown\n";
      exit(1);
  }
}
