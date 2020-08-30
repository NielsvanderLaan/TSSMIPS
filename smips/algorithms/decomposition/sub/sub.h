#ifndef SUB_H
#define SUB_H

#include "gurobi_c++.h"

#include "../../../problem_data/problem.h"

using namespace std;

class Sub
{
  public:
    GRBModel d_model;
    Problem &d_problem;
    size_t d_m2, d_n2;
    vector<GRBConstr> d_constrs;
    vector<GRBVar> d_vars;
    vector<double> &d_q;
    
    Sub(GRBEnv &env, Problem &problem);
    Sub(const Sub &other);
    ~Sub();

    void update(double *rhs, size_t s);
    vector<double> compute_slope(size_t, double *x);
    
    struct GomInfo
    {
      double *lambda;
      int *vBasis;
      int *cBasis;
    };
    
    struct Multipliers
    {
      double *lambda;
      double obj;
    };
    
    
    Multipliers solve();
    GomInfo solve2();
};

#endif
