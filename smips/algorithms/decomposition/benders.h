#ifndef BENDERS_H
#define BENDERS_H

#include <algorithm>
#include <chrono>
#include <string>

#include "../../problem_data/problem.h"
#include "master/master.h"
#include "sub/sub.h"
#include "lagrangian/lagrangian.h"
#include "gomory/gomory.h"
#include "ald/ald.h"
#include "pslp/pslp.h"
#include "aggregator/aggregator.h"
#include "cut/benderscut.h"
#include <memory>

using namespace std;

class Benders
{
  public:
    Problem &d_problem;  // contains problem data
    GRBEnv &d_env;
    size_t d_n1, d_p1, d_m2, d_n2, d_S;
    Master d_master;    // master problem
    Sub d_sub;          // sub-problem
    Lagrangian d_lr;    // lagrangian relaxation   
    Gomory d_gomory;    // Gomory relaxation  
    Ald d_ald;          // For deriving ALD cuts  
    Pslp d_pslp;        // For deriving (strong) ZK cuts  
    Aggregator d_agg;
    
    vector<double> d_lb, d_ub;
    vector<vector<vector<double>>> d_visited;  // for each scenario, we store the basis matrices that we have visited (encoded by vBasis, cBasis)
    vector<vector<double>> d_objectives;       // for each visited basis matrix, we store the corresponding gomory objective value   
    
    double *d_xvals;
    double *d_incumbent;    // keeps track of best solution
    
    void computeTx(double *x, double *Tx);     // computes Tx (rba)

    Benders(GRBEnv &env, GRBenv *c_env, Problem &problem, bool zk_safe = true);    // initializes d_master and d_sub with both arguments  
    Benders(const Benders &other);
                                               
    ~Benders();
    
    double get_lb();
    void update_bounds(size_t var, double val, bool lower);

    BendersCut lpCut(double *x); 
    BendersCut sb_cut(double *x);     
    void ald_cut(double *x, double *beta, double &tau, double &gamma, size_t maxRounds); // RBA
    BendersCut lbdaCut(double *x, double *alpha);
    double compute_gomory(size_t s, int *vBasis, int *cBasis, double *ws, double *alpha);
    
    struct Bounds { double d_LB; double d_UB; };
    
    double lpSolve(double tol = 1e-4);                                               // L_shaped 
    void strong_benders(double tol = 1e-4);                                          // uses strengthened L-shaped cuts
    void lbda(double *alpha, double gomoryTimeLimit = 1e6, double tol = 1e-4);       // LBDA(alpha)  
    void ald_solve(double tol = 1e-4, size_t maxRounds = 25);    
    double zk_solve(double tol = 1e-4, size_t maxRounds = 25);    
    Bounds hybrid_solve(double global_UB = GRB_INFINITY, bool affine = false, bool lp_cuts= false, double tol = 1e-4);
    

    //bool add_cut(double *beta, double gamma, double kappa, double *x , double theta, double tol); 
    bool add_cut(BendersCut &cut, Master::Solution sol, double tol);
              // adds the cut kappa theta - beta^T x>= gamma  and returns false if cut was added
              // add_cut() updates the master problem, but also the cglp objects via pslp
              
    bool round_of_cuts(Master::Solution sol, double tol); 
};

#endif












