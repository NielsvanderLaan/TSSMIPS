#ifndef TSSMIPS_RUN_H
#define TSSMIPS_RUN_H

#include <iostream>
#include <chrono>
#include <vector>

#include <omp.h>
#include <cstdio>
#include "assert.h"
#include "gurobi_c++.h"
#include "gurobi_c.h"

#include "../smips/problem_data/problem.h"
#include "../smips/algorithms/deqform/deqform.h"
#include "../smips/algorithms/decomposition/benders.h"
#include "../smips/algorithms/trees/tree.h"

void instance(Problem &problem, int argc, char *argv[]);
bool find(int argc, char *argv[], string cmp);
bool use_rcuts(int argc, char *argv[]);
bool use_fenchel(int argc, char *argv[]);
bool solve_DEF(int argc, char *argv[]);
bool solve_root(int argc, char *argv[]);
bool solve_tree(int argc, char *argv[]);
int get_max_rounds(int argc, char *argv[]);
int nthreads(int argc, char* argv[]);
double get_time_limit(int argc, char *argv[]);
void setLogger(string filename);


void details(vector<Type> types, int max_rounds, bool rcuts, bool fenchel, double time_limit, int thread_count);

#endif //TSSMIPS_RUN_H
