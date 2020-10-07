#ifndef BENDERS_H
#define BENDERS_H

#include <algorithm>
#include <chrono>
#include <string>
#include <list>
#include <numeric>

#include "type.h"
#include "../../../debug.h"
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
    size_t d_n1, d_p1, d_n2, d_m2, d_S;
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
    double d_UB;
    
    void computeTx(double *x, double *Tx);     // computes Tx (rba)

    Benders(GRBEnv &env, GRBenv *c_env, Problem &problem, bool zk_safe = true);    // initializes d_master and d_sub with both arguments  
    Benders(const Benders &other);
                                               
    ~Benders();
    
    double get_lb();
    void update_bounds(size_t var, double val, bool lower);

    BendersCut lpCut(double *x); 
    BendersCut sb_cut(double *x);
    BendersCut lr_cut(double *x, vector<double> &vx);
    void ald_cut(double *x, double *beta, double &tau, double &gamma, size_t maxRounds); // RBA
    BendersCut lbdaCut(double *x, double *alpha);
    double compute_gomory(size_t s, int *vBasis, int *cBasis, double *ws, double *alpha);
    
    struct Bounds { double d_LB; double d_UB; bool branch; };
    
    double lpSolve(double tol = 1e-4);                                               // L_shaped
    void lbda(double *alpha, double gomoryTimeLimit = 1e6, double tol = 1e-4);       // LBDA(alpha)  
    void ald_solve(double tol = 1e-4, size_t maxRounds = 25);
    Bounds hybrid_solve(vector<Type> types, bool force_int, size_t max_rounds = 25,
                        double upper_bound = GRB_INFINITY, double tol = 1e-4, double time_limit = 1e100, bool rcuts = true);

    //bool add_cut(double *beta, double gamma, double kappa, double *x , double theta, double tol); 
    bool add_cut(BendersCut &cut, Master::Solution sol, double tol);
              // adds the cut kappa theta - beta^T x>= gamma  and returns false if cut was added
              // add_cut() updates the master problem, but also the cglp objects via pslp
    void reverse_cut(double UB);
    void update(double UB, bool rcuts = true);
              
    size_t round_of_cuts(Master::Solution sol, double tol);

    BendersCut compute_cut(Type type, Master::Solution &sol, bool int_feas, vector<double> &vx, double tol, double *alpha = nullptr);

};

static double avg(double time, size_t n) {return n == 0 ? 0.0 : time / n;}

#endif












