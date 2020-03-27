#include <iostream>
#include <chrono>
#include <vector>

#include "gurobi_c++.h"
#include "gurobi_c.h"

#include "smips/problem_data/problem.h"
#include "smips/algorithms/deqform/deqform.h"
#include "smips/algorithms/decomposition/benders.h"
#include "smips/algorithms/trees/tree.h"

using namespace std;

int main(int argc, char *argv[])
{    
  Data rand(46511);

  GRBEnv env;  
  env.set(GRB_IntParam_OutputFlag, 0); 
  env.set(GRB_IntParam_Threads, 1);
  
  GRBenv *c_env;
  GRBloadenv(&c_env, NULL);
  GRBsetintparam(c_env, "OutputFlag", 0);
  GRBsetintparam(c_env, "Threads", 1);

  {  
    size_t n1, p1, m1, n2, p2, m2, S;            // input size
    
    n1 = 10; p1 = 10; m1 = 4; n2 = 12; p2 = 12; m2 = 6; S = 1000;
                                                 // parameter bounds (uniform distribution)  
    size_t A_low, A_high, T_low, T_high, W_low, W_high, c_low, c_high, b_low, b_high, q_low, q_high;
    A_low = 1; A_high = 4; T_low = 1; T_high = 3; W_low = 1; W_high = 2; 
    c_low = 1; c_high = 3; b_low = 10, b_high = 15; q_low = 25; q_high = 35;
    
                                                 // create problem
    Problem problem(n1, p1, m1, n2, p2, m2, S, rand, env, m1, 0, 0, m2); 
    
    problem.randomInstance(A_low, A_high, T_low, T_high, W_low, W_high, c_low, c_high, b_low, b_high, q_low, q_high);
    problem.set_omega_gaus(25.0, 2.0);   
    problem.enforce_ccr(50.0);
    
    vector<double> l1(n1, 0.0); vector<double> u1(n1, 5.0); vector<double> l2(n2, 0.0); vector<double> u2(n2, 20.0); 
    problem.set_bounds(l1, u1, l2, u2);

    double *x;

    Tree tree(env, c_env, problem);

    auto t1 = chrono::high_resolution_clock::now();
    vector<double> x_bab = tree.bab(false, 1e-2);
    auto t2 = chrono::high_resolution_clock::now();
    for_each(x_bab.begin(), x_bab.end(), [](double val){cout << val << ' ';});
    cout << "\ncx + Q(x) = " << problem.evaluate(x_bab.data()) << '\n';
    cout << "computation time: " << chrono::duration_cast<chrono::milliseconds>(t2 - t1).count() / 1000.0 << '\n';

    /*
    Benders ben(env, c_env, problem); 
    cout << "-------------L-shaped------------\n"; \
       
    cout << "lp lower bound: " << ben.lpSolve() << '\n';    
    x = ben.d_xvals; 
    cout << "x: ";  
    for (size_t var = 0; var != n1; ++var)
      cout << x[var] << ' ';
    cout << '\n';
    cout << "cx + Q(x) = " << problem.evaluate(x) << '\n'; 
    

    auto t1 = chrono::high_resolution_clock::now();
    ben.hybrid_solve(1e-2);
    auto t2 = chrono::high_resolution_clock::now();

    
    //cout << "zk_solve() lb: " << ben.zk_solve() << '\n';
    */
    
    
    /*
    x = ben.d_xvals; 
    cout << "x: ";  
    for (size_t var = 0; var != n1; ++var)
      cout << x[var] << ' ';
    cout << '\n';
    cout << "cx + Q(x) = " << problem.evaluate(x) << '\n'; 
    */
    
    //cout << "computation time: " << chrono::duration_cast<chrono::milliseconds>(t2 - t1).count() / 1000.0 << '\n';  


  
 
 
 
    /*
    cout << "-------------ZK cuts------------\n";
    cout << "zk lower bound = " << ben.zk_solve(0.1) << '\n';
    
    x = ben.d_xvals;
    cout << "x: ";  
    for (size_t var = 0; var != n1; ++var)
      cout << x[var] << ' ';
    cout << '\n';
    cout << "cx + Q(x) = " << problem.evaluate(x) << '\n';
    */
    

    cout << "-------------Solving DEF------------\n";

    DeqForm DEF(env, problem);   
    DEF.d_model.set(GRB_IntParam_OutputFlag, 1);
    DEF.solve(1500.0);
          // printing DEF solution and objective
    x = DEF.d_xVals;
    cout << "x: ";
    for (size_t var = 0; var != n1; ++var)
      cout << x[var] << ' ';  
    cout << '\n'; 
    cout <<  "cx + Q(x) = " << problem.evaluate(x) << '\n';

  } 

  GRBfreeenv(c_env);
}












