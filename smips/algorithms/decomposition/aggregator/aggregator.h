#ifndef AGG_H
#define AGG_H

#include "cgmip/cgmip.h"
#include "../master/master.h"

class Aggregator
{
  public:
    size_t d_n1;
    vector<double> d_probs;
    vector<CGMip> d_cgmips;
        // for computing v_w(x) and Q(x)
    GRBModel d_vw;        // used to evaluate the value function
    vector<vector<double>> &d_omega;
    vector<vector<double>> &d_Tmat;
    vector<vector<double>> &d_q_omega;
    vector<vector<vector<double>>> d_W_omega;
    bool d_fix_rec;
    
    Aggregator(GRBEnv &env, Problem &problem);
    
    void add_rows(BendersCut &cut);
    double compute_vwx(double *x, size_t s);
    BendersCut strong_cut(Master::Solution sol, double &Qx, bool affine, double tol, double rho_tol = 1e-4);
    
    
    
    void update_bounds(size_t var, double val, bool lower);
};

#endif