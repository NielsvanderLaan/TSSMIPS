#include "sub.h"

Sub::Sub(GRBEnv &env, Problem &problem, size_t scenario)
:
d_model(env),
d_q(problem.d_q)
{
  size_t n2 = problem.d_n2;
  size_t m2 = problem.d_m2;
  
  size_t ss_leq, ss_geq;
  ss_leq = problem.d_ss_leq; ss_geq = problem.d_ss_geq;
  
  d_n2 = n2;
  d_m2 = m2;

       // variable types    
  char vTypes[n2];
  fill(vTypes, vTypes + n2, GRB_CONTINUOUS);    
      // cost vector
  double *q;
  if (scenario == -1)
    q = problem.d_q.data();        // transform cost vector and omega to c-style array
  else
    q = problem.d_q_omega[scenario].data();
      // add variables

  GRBVar* var_ptr = d_model.addVars(problem.d_l2.data(), problem.d_u2.data(), q, vTypes, NULL, n2);
  for (size_t var = 0; var != n2; ++var)
    d_vars.push_back(var_ptr[var]);
  delete[] var_ptr;
      // constraint senses
  char senses[m2];
  fill(senses,                   senses + ss_leq,          GRB_LESS_EQUAL);
  fill(senses + ss_leq,          senses + ss_leq + ss_geq, GRB_GREATER_EQUAL);
  fill(senses + ss_leq + ss_geq, senses + m2,              GRB_EQUAL);
  
  
      // constraint rhs
  double rhs[m2];
  fill_n(rhs, m2, 0.0);
  

      // constraint lhs
  vector<vector<double>> &Wmat = scenario == -1 ? problem.d_Wmat : problem.d_W_omega[scenario];
  GRBLinExpr Wy[m2];
  for (size_t conIdx = 0; conIdx != m2; ++conIdx)
  {
    double *row = Wmat[conIdx].data();      
    Wy[conIdx].addTerms(row, d_vars.data(), n2);
  }
      // add constraints
  GRBConstr *constr_ptr = d_model.addConstrs(Wy, senses, rhs, NULL, m2);
  for (size_t con = 0; con != m2; ++con)
    d_constrs.push_back(constr_ptr[con]);
  delete[] constr_ptr;
  d_model.update();
}