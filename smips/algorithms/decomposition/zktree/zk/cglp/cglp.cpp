#include "cglp.h"

Cglp::Cglp(Problem &problem, GRBEnv &env, size_t scenario)
:  
  d_model(env),
  d_n1(problem.d_n1),
  d_p1(problem.d_p1),
  d_m1(problem.d_m1),
  d_n2(problem.d_n2), 
  d_p2(problem.d_p2), 
  d_m2(problem.d_m2),
  d_l1_mults(d_n1, -1),
  d_u1_mults(d_n1, -1),
  d_l2_mults(d_n2, -1),
  d_u2_mults(d_n2, -1)
{  
  //d_model.set(GRB_IntParam_OutputFlag, 1);
  d_model.set(GRB_IntParam_ScaleFlag, 0);
      // adding cut coefficients as decision variables (normalization imposed here)
  for (size_t var = 0; var != d_n1; ++var)
    d_Trow.push_back(d_model.addVar(-1, 1, 0.0, GRB_CONTINUOUS));    
  for (size_t var = 0; var != d_n2; ++var)
    d_Wrow.push_back(d_model.addVar(-1, 1, 0.0, GRB_CONTINUOUS));

  d_r = d_model.addVar(0, 1, 0, GRB_CONTINUOUS);
  d_h = d_model.addVar(-1, 1, -1, GRB_CONTINUOUS);
  
      // adding the multipliers
  size_t fs_eq = d_m1 - problem.d_fs_leq - problem.d_fs_geq;
  size_t ss_eq = d_m2 - problem.d_ss_leq - problem.d_ss_geq;
  
  d_lambda1.push_back(d_model.addVar(-GRB_INFINITY, 0, 0, GRB_CONTINUOUS));    // corresponding to the disjunction
  d_lambda2.push_back(d_model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS));
  
  for (size_t con = 0; con != problem.d_fs_leq; ++con)  // first-stage <= constraints
  {  
    d_lambda1.push_back(d_model.addVar(-GRB_INFINITY, 0, 0, GRB_CONTINUOUS));
    d_lambda2.push_back(d_model.addVar(-GRB_INFINITY, 0, 0, GRB_CONTINUOUS));
  }
  for (size_t con = 0; con != problem.d_fs_geq; ++con)  // first-stage >= constraints
  {
    d_lambda1.push_back(d_model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS));
    d_lambda2.push_back(d_model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS));
  }
  for (size_t con = 0; con != fs_eq; ++con)  // first-stage == constraints
  {
    d_lambda1.push_back(d_model.addVar(-GRB_INFINITY, GRB_INFINITY, 0, GRB_CONTINUOUS));
    d_lambda2.push_back(d_model.addVar(-GRB_INFINITY, GRB_INFINITY, 0, GRB_CONTINUOUS));
  }
  for (size_t con = 0; con != problem.d_ss_leq; ++con)  // second-stage <= constraints
  {
    d_lambda1.push_back(d_model.addVar(-GRB_INFINITY, 0, 0, GRB_CONTINUOUS));
    d_lambda2.push_back(d_model.addVar(-GRB_INFINITY, 0, 0, GRB_CONTINUOUS));
  }
  for (size_t con = 0; con != problem.d_ss_geq; ++con)  // second-stage >= constraints
  {
    d_lambda1.push_back(d_model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS));
    d_lambda2.push_back(d_model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS));
  }
  for (size_t con = 0; con != ss_eq; ++con)  // second-stage == constraints
  {
    d_lambda1.push_back(d_model.addVar(-GRB_INFINITY, GRB_INFINITY, 0, GRB_CONTINUOUS));
    d_lambda2.push_back(d_model.addVar(-GRB_INFINITY, GRB_INFINITY, 0, GRB_CONTINUOUS));
  }
  
  d_nMults = d_lambda1.size();  

      // adding the constraints
  for (size_t var = 0; var != d_n1; ++var)          // looping over x variables
  {
    vector<double> coeffs(d_nMults);
    for (size_t row = 0; row != d_m1; ++row)
      coeffs[1 + row] = problem.d_Amat[row][var];
    for (size_t row = 0; row != d_m2; ++row)
      coeffs[1 + d_m1 + row] = problem.d_Tmat[row][var];
    
    GRBLinExpr lhs1;
    lhs1.addTerms(coeffs.data(), d_lambda1.data(), d_nMults); 
    d_constrs1.push_back(d_model.addConstr(lhs1, GRB_LESS_EQUAL, d_Trow[var]));
    
    GRBLinExpr lhs2;
    lhs2.addTerms(coeffs.data(), d_lambda2.data(), d_nMults); 
    d_constrs2.push_back(d_model.addConstr(lhs2, GRB_LESS_EQUAL, d_Trow[var]));
  }
  
  d_constrs1.push_back(d_model.addConstr(0.0, GRB_LESS_EQUAL, d_r));    // c.t. theta
  d_constrs2.push_back(d_model.addConstr(0.0, GRB_LESS_EQUAL, d_r));
  
  for (size_t var = 0; var != d_n2; ++var)             // looping over y-variables
  {
    vector<double> coeffs(d_nMults);
    for (size_t row = 0; row != d_m2; ++row)
      coeffs[1 + d_m1 + row] = problem.d_Wmat[row][var];
  
    GRBLinExpr lhs1;
    lhs1.addTerms(coeffs.data(), d_lambda1.data(), d_nMults); 
    d_constrs1.push_back(d_model.addConstr(lhs1, GRB_LESS_EQUAL, d_Wrow[var]));
    
    GRBLinExpr lhs2;
    lhs2.addTerms(coeffs.data(), d_lambda2.data(), d_nMults); 
    d_constrs2.push_back(d_model.addConstr(lhs2, GRB_LESS_EQUAL, d_Wrow[var]));
  }
  
  double coeffs[d_nMults];                            // c.t. rhs
  coeffs[0] = 0;
  copy_n(problem.d_b.data(),  d_m1, coeffs + 1);
  copy_n(problem.d_omega[scenario].data(), d_m2, coeffs + 1 + d_m1);
  
  GRBLinExpr lhs1;
  lhs1.addTerms(coeffs, d_lambda1.data(), d_nMults); 
  d_constrs1.push_back(d_model.addConstr(lhs1, GRB_GREATER_EQUAL, d_h));
  
  GRBLinExpr lhs2;
  lhs2.addTerms(coeffs, d_lambda2.data(), d_nMults); 
  d_constrs2.push_back(d_model.addConstr(lhs2, GRB_GREATER_EQUAL, d_h));
 
  for (size_t var = 0; var != d_n1; ++var)      // first-stage bounds
  {
    GRBConstr constrs1[2] = {d_constrs1[var], d_constrs1[d_n1 + d_n2 + 1]};
    GRBConstr constrs2[2] = {d_constrs2[var], d_constrs2[d_n1 + d_n2 + 1]};
      double lb = problem.d_l1[var];
    if (lb > 0)
    {
      double coeffs[2] = {1, lb};
      d_lambda1.push_back(d_model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS, 2, constrs1, coeffs));
      d_lambda2.push_back(d_model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS, 2, constrs2, coeffs));
      d_l1_mults[var] = d_nMults;
      ++d_nMults;
    }
    
    double ub = problem.d_u1[var];
    if (ub < 1e10)
    {
      double coeffs[2] = {1, ub};
      d_lambda1.push_back(d_model.addVar(-GRB_INFINITY, 0, 0, GRB_CONTINUOUS, 2, constrs1, coeffs));
      d_lambda2.push_back(d_model.addVar(-GRB_INFINITY, 0, 0, GRB_CONTINUOUS, 2, constrs2, coeffs));
      d_u1_mults[var] = d_nMults;
      ++d_nMults;
    }
  }

  for (size_t var = 0; var != d_n2; ++var)      // second-stage bounds
  {
    GRBConstr constrs1[2] = {d_constrs1[d_n1 + 1 + var], d_constrs1[d_n1 + d_n2 + 1]};
    GRBConstr constrs2[2] = {d_constrs2[d_n1 + 1 + var], d_constrs2[d_n1 + d_n2 + 1]};
    double lb = problem.d_l2[var];
    if (lb > 0)
    {
      double coeffs[2] = {1, lb};
      d_lambda1.push_back(d_model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS, 2, constrs1, coeffs));
      d_lambda2.push_back(d_model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS, 2, constrs2, coeffs));
      d_l2_mults[var] = d_nMults;
      ++d_nMults;
    }
    
    double ub = problem.d_u2[var];
    if (ub < 1e10)
    {
      double coeffs[2] = {1, ub};
      d_lambda1.push_back(d_model.addVar(-GRB_INFINITY, 0, 0, GRB_CONTINUOUS, 2, constrs1, coeffs));
      d_lambda2.push_back(d_model.addVar(-GRB_INFINITY, 0, 0, GRB_CONTINUOUS, 2, constrs2, coeffs));
      d_u2_mults[var] = d_nMults;
      ++d_nMults;
    }
  } 

  
  d_model.update();  
}


















