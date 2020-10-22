#ifndef AGG_H
#define AGG_H

#include "cgmip/cgmip.h"
#include "../zktree/zktree.h"
#include "../master/master.h"

class Aggregator
{
  public:
    size_t d_n1;
    vector<double> d_probs;
    vector<CGMip> d_cgmips;
    vector<ZkTree> d_trees;

    // for computing v_w(x) and Q(x)
    GRBModel d_vw;        // used to evaluate the value function
    Problem &d_problem;
    bool d_fix_rec;
    Aggregator(GRBEnv &env, GRBenv *c_env, Problem &problem);
    
    void add_rows(BendersCut &cut);
    vector<double> compute_vwx(double *x);
    BendersCut strong_cut(Master::Solution sol, vector<double> &vx, bool affine, double tol, bool int_feas = true, double rho_tol = 1e-4);
    BendersCut bac_cut(Master::Solution sol, Master &mp, double tol, size_t maxRounds = 25, double rho_tol = 1e-4);

    void update_bounds(size_t var, double val, bool lower);
    void reverse_cut(double UB);
};

#endif