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
    vector<vector<double>> &d_omega;
    vector<vector<double>> &d_Tmat;
    vector<vector<double>> &d_q_omega;
    vector<vector<vector<double>>> &d_W_omega;
    bool d_fix_rec;
    
    Aggregator(GRBEnv &env, GRBenv *c_env, Problem &problem);
    
    void add_rows(BendersCut &cut);
    vector<double> compute_vwx(double *x);
    BendersCut strong_cut(Master::Solution sol, vector<double> &vx, bool affine, double tol, double rho_tol = 1e-4);
    BendersCut bac_cut(Master::Solution sol, Master &mp, double tol, size_t maxRounds = 25, double rho_tol = 1e-4);

    void update_bounds(size_t var, double val, bool lower);
};

#endif