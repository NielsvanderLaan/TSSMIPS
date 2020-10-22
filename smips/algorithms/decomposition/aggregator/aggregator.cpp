#include "aggregator.h"

Aggregator::Aggregator(GRBEnv &env, GRBenv *c_env, Problem &problem)
:  
  d_n1(problem.d_n1),
  d_probs(problem.d_probs),
  d_vw(env),
  d_problem(problem),
  d_fix_rec(problem.d_fix_rec)
{
  d_cgmips.reserve(problem.d_S);
  d_trees.reserve(problem.d_S);

  for (size_t s = 0; s != problem.d_S; ++s)
  {
    CGMip cgmip{ env, problem, s };
    d_cgmips.push_back(cgmip);

    ZkTree tree{ c_env, env, problem, s };
    d_trees.push_back(tree);
  }

      // Initializing d_vw
  size_t n2 = problem.d_n2;
  size_t p2 = problem.d_p2;    
  size_t m2 = problem.d_m2; 
      
  vector<char> y_types(n2, GRB_CONTINUOUS);
  fill_n(y_types.begin(), p2, GRB_INTEGER);        
  GRBVar *vw_vars = d_vw.addVars(problem.d_l2.data(), problem.d_u2.data(), problem.d_q.data(), y_types.data(), NULL, n2);
  GRBLinExpr Wy[m2];
  for (size_t con = 0; con != m2; ++con)
    Wy[con].addTerms(problem.d_Wmat[con].data(), vw_vars, n2);
  vector<double> rhs(m2);
  vector<char> senses(m2, GRB_EQUAL);
  fill_n(senses.begin(), problem.d_ss_leq, GRB_LESS_EQUAL);
  fill_n(senses.begin() + problem.d_ss_leq, problem.d_ss_geq, GRB_GREATER_EQUAL);
  delete[] d_vw.addConstrs(Wy, senses.data(), rhs.data(), NULL, m2);
  delete[] vw_vars;
  
  d_vw.update();  
}