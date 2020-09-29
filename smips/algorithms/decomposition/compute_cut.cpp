#include"benders.h"

BendersCut Benders::compute_cut(Type type, Master::Solution &sol, bool int_feas, vector<double> &vx, double tol, double *alpha)
{
  switch (type)
  {
    case LP:        return lpCut(sol.xVals.data());
    case SB:        return sb_cut(sol.xVals.data());
    case LR:        return d_agg.strong_cut(sol, vx, true, tol, int_feas);
    case LR_LAP:    return d_pslp.best_zk_cut(sol,d_master, true, true);
    case SC_ZK:     return d_pslp.best_zk_cut(sol,d_master);
    case SC_LAP:    return d_pslp.best_zk_cut(sol,d_master, true);
    case SC_BAB:    return d_agg.bac_cut(sol, d_master, tol);
    case SC_RG:     return d_agg.strong_cut(sol, vx, false, tol, int_feas);
    case LBDA:
      if (not alpha)
      {
        cerr << "alpha = nullptr\n";
        exit(1);
      }
      return lbdaCut(sol.xVals.data(), alpha);
    default:
      cerr << "cut type unknown\n";
      exit(1);
  }
}
