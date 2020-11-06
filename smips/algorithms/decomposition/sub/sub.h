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
    int d_s;
    size_t d_m2, d_n2;
    vector<GRBConstr> d_constrs;
    vector<GRBVar> d_vars;

    Sub(GRBEnv &env, Problem &problem, size_t s);
    Sub(const Sub &other);

    void update(vector<double> &x);
    vector<double> compute_slope(vector<double> &x);
    
    struct GomInfo
    {
      vector<double> lambda;
      vector<int> vbasis;
      vector<int> cbasis;
    };
    
    struct Multipliers
    {
      vector<double> lambda;
      double obj;
    };
    
    
    Multipliers solve();
    GomInfo solve2();
};


#endif
