#ifndef TSSMIPS_RUN_H
#define TSSMIPS_RUN_H

#include <iostream>
#include <chrono>
#include <vector>

#include "gurobi_c++.h"
#include "gurobi_c.h"

#include "../smips/problem_data/problem.h"
#include "../smips/algorithms/deqform/deqform.h"
#include "../smips/algorithms/decomposition/benders.h"
#include "../smips/algorithms/trees/tree.h"

void solve_ssv(Data &rand, GRBEnv &env, GRBenv *c_env, vector<Type> types, bool rcuts);
void lbda_scheme(Problem &problem, GRBEnv &env, GRBenv *c_env, size_t nIter);
void solve_dcap(vector<Type> types, Data &rand, GRBEnv &env, GRBenv *c_env);
void solve_sizes(vector<Type> types, Data &rand, GRBEnv &env, GRBenv *c_env);
void solve_ri(Data &rand, GRBEnv &env, GRBenv *c_env);
void solve_caroe(Data &rand, GRBEnv &env, GRBenv *c_env);

#endif //TSSMIPS_RUN_H
