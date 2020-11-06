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
#include "gomory/gomory.h"
#include "aggregator/aggregator.h"
#include "cut/benderscut.h"

using namespace std;

class Benders
{
  public:
    Problem &d_problem;  // contains problem data
    GRBEnv &d_env;
    size_t d_n1, d_p1, d_n2, d_m2, d_S;
    Master d_master;    // master problem
    Gomory d_gomory;    // Gomory relaxation
    Aggregator d_agg;
    
    vector<double> d_lb, d_ub;
    vector<vector<vector<double>>> d_visited;  // for each scenario, we store the basis matrices that we have visited (encoded by vBasis, cBasis)
    vector<vector<double>> d_objectives;       // for each visited basis matrix, we store the corresponding gomory objective value   
    
    double *d_xvals;
    double *d_incumbent;    // keeps track of best solution
    double d_UB;


    Benders(GRBEnv &env, GRBenv *c_env, Problem &problem, bool zk_safe = true);    // initializes d_master and d_sub with both arguments  
    Benders(const Benders &other);
                                               
    ~Benders();
    
    double get_lb();
    void update_bounds(size_t var, double val, bool lower);

    BendersCut lbdaCut(vector<double> &x, vector<double> &alpha);
    double compute_gomory(size_t s, vector<int> &vbasis, vector<int> &cbasis, vector<double> &ws, vector<double> &alpha);
    
    struct Bounds { double d_LB; double d_UB; bool branch; };
    
    double lpSolve(double tol = 1e-4);                                               // L_shaped
    void lbda(vector<double> &alpha, double gomoryTimeLimit = 1e6, double tol = 1e-4);       // LBDA(alpha)

    Bounds hybrid_solve(vector<Type> types, bool force_int, int max_rounds = 25,
                        double upper_bound = GRB_INFINITY, double tol = 1e-4, double time_limit = 1e100,
                        bool rcuts = true, bool fenchel = true);

    //bool add_cut(double *beta, double gamma, double kappa, double *x , double theta, double tol); 
    bool add_cut(BendersCut &cut, Master::Solution sol, double tol);
              // adds the cut kappa theta - beta^T x>= gamma  and returns false if cut was added
              // add_cut() updates the master problem, but also the cglp objects via pslp
    void reverse_cut(double UB);
    void update(double UB, bool rcuts);
              
    size_t round_of_cuts(Master::Solution sol, double tol);

    BendersCut compute_cut(Type type, Master::Solution &sol, bool int_feas, vector<double> &vx, double tol, vector<double> alpha = vector<double>(0));
};

static double avg(double time, size_t n) {return n == 0 ? 0.0 : time / n;}

#endif












