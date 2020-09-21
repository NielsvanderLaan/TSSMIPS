#ifndef PROBLEM_H
#define PROBLEM_H

#include <vector>
#include <algorithm>

#include <fstream>    // for reading data

#include "gurobi_c++.h"
#include "data.h"

class Problem
{
  public:
    size_t d_n1, d_p1, d_m1, d_n2, d_p2, d_m2, d_S;  // size parameters
    size_t d_fs_leq, d_fs_geq, d_ss_leq, d_ss_geq;       // number of >= and <= constraints in the first and second stage
    bool d_fix_rec;
    vector<double> d_l1, d_u1, d_l2, d_u2;
    vector<double> d_c, d_b, d_q;
    vector<vector<double>> d_Amat, d_Tmat, d_Wmat, d_omega;   // rows of d_omega correspond to scenarios
    vector<vector<double>> d_q_omega;
    vector<vector<vector<double>>> d_W_omega;
    vector<double> d_probs;
    double d_L; // lb of v
    
    Data &d_gen; // used to generate (random) data
    GRBEnv &d_env;
    
    Problem(Data &generator, GRBEnv &env);
    Problem(size_t n1, size_t p1, size_t m1, 
            size_t n2, size_t p2, size_t m2, size_t S,
            Data &generator, GRBEnv &env,
            size_t fs_leq = 0, size_t fs_geq = 0,
            size_t ss_leq = 0, size_t ss_geq = 0); // initializes size parameters, d_gen, and d_sub
    Problem(const Problem &other) = delete;
    ~Problem();
          
    void randomInstance(double sigma = 0.1, double mu = 10.0,
                        int A_low = 1, int A_high = 6,
                        int T_low = 1, int T_high = 6,
                        int W_low = 1, int W_high = 6,
                        int c_low = 1, int c_high = 5,
                        int b_low = 1, int b_high = 5,
                        int q_low = 5, int q_high = 10);   // initializes A, T, W, b, c, q (>= constraints)
                        
    void enforce_ccr(double penalty);
  
    void set_omega_gaus(double mean, double sd);    // initializes d_omega  
    void set_bounds(vector<double> &l1, vector<double> &u1, 
                    vector<double> &l2, vector<double> &u2);
                    
    void sizes(size_t S);
    void ssv95(size_t S, bool fs_continuous, bool ss_binary, bool standard_T = true);
    void sslp(size_t nServers, size_t nClients, size_t S);
    void dcap(size_t nResources, size_t nClients, size_t nPeriods, size_t S, bool fs_cont = false, double fs_ub = 1.0);
    void classic_ri();
    void caroe(size_t S);
    void caroe_LD(size_t S);

    
    double evaluate(double *x);  // evaluates cx + Q(x) (does not check feasibility)      
    

    GRBModel* init_sub();
};

#endif