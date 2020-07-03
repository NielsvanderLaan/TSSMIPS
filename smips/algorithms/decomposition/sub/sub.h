#ifndef SUB_H
#define SUB_H

#include "gurobi_c++.h"

#include "../../../problem_data/problem.h"

using namespace std;

class Sub
{
  public:
    GRBModel d_model;
    size_t d_m2, d_n2;
    vector<GRBConstr> d_constrs;
    vector<GRBVar> d_vars;
    vector<double> &d_q;
    
    Sub(GRBEnv &env, Problem &problem, size_t scenario = -1);
    Sub(const Sub &other);
    ~Sub();

    void update(double *rhs);
    
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
