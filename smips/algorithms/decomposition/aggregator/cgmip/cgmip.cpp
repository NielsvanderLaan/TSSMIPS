#include "cgmip.h"

CGMip::CGMip(GRBEnv &env, Problem &problem, size_t s)
:
  d_mp(env),
  d_sub(env),
  d_beta(problem.d_n1),
  d_xVars(problem.d_n1),
  d_yVars(problem.d_n2)
{ 
  size_t n1 = problem.d_n1;
  size_t p1 = problem.d_p1;
  size_t m1 = problem.d_m1;
  size_t n2 = problem.d_n2;
  size_t p2 = problem.d_p2;
  size_t m2 = problem.d_m2;
 
  // Initializing MP
  vector<double> lb(n1, -1e20);
   
  d_alpha = d_mp.addVar(-1e20, 1e20, -1.0, GRB_CONTINUOUS);
  GRBVar *beta = d_mp.addVars(lb.data(), NULL, NULL, NULL, NULL, n1);
  d_tau = d_mp.addVar(0.0, 1e20, 0.0, GRB_CONTINUOUS);

  // Initializing SP
    // adding xvars
  vector<char> x_types(n1, GRB_CONTINUOUS);
  fill_n(x_types.begin(), p1, GRB_INTEGER);
  GRBVar *xVars = d_sub.addVars(problem.d_l1.data(), problem.d_u1.data(), NULL, x_types.data(), NULL, n1);
    // adding theta and eta
  d_theta = d_sub.addVar(0, 1e20, 0, GRB_CONTINUOUS);
  d_eta = d_sub.addVar(0, 1e20, 1.0, GRB_CONTINUOUS);
    // adding yvars  
  vector<char> y_types(n2, GRB_CONTINUOUS);
  fill_n(y_types.begin(), p2, GRB_INTEGER);
  GRBVar *yVars = d_sub.addVars(problem.d_l2.data(), problem.d_u2.data(), NULL, y_types.data(), NULL, n2);
    // adding Ax ~ b
  GRBLinExpr Ax[m1];
  for (size_t con = 0; con != m1; ++con)
    Ax[con].addTerms(problem.d_Amat[con].data(), xVars, n1);
  vector<char> senses1(m1, GRB_EQUAL);
  fill_n(senses1.begin(), problem.d_fs_leq, GRB_LESS_EQUAL);
  fill_n(senses1.begin() + problem.d_fs_leq, problem.d_fs_geq, GRB_GREATER_EQUAL);
  delete[] d_sub.addConstrs(Ax, senses1.data(), problem.d_b.data(), NULL, m1); 
    // adding Wy + Tx ~ omega
  GRBLinExpr WyTx[m2];
  for (size_t con = 0; con != m2; ++con)
  {
    WyTx[con].addTerms(problem.d_Wmat[con].data(), yVars, n2);
    WyTx[con].addTerms(problem.d_Tmat[con].data(), xVars, n1);
  }
  vector<char> senses2(m2, GRB_EQUAL);
  fill_n(senses2.begin(), problem.d_ss_leq, GRB_LESS_EQUAL);
  fill_n(senses2.begin() + problem.d_ss_leq, problem.d_ss_geq, GRB_GREATER_EQUAL);
  delete[] d_sub.addConstrs(WyTx, senses2.data(), problem.d_omega[s].data(), NULL, m2);
    // adding eta >= q^T y
  GRBLinExpr qy;
  qy.addTerms(problem.d_q.data(), yVars, n2);
  d_sub.addConstr(d_eta, GRB_GREATER_EQUAL, qy);
  

     // cleaning up
  copy_n(beta, n1, d_beta.begin());
  copy_n(xVars, n1, d_xVars.begin());
  copy_n(yVars, n2, d_yVars.begin());
  delete[] beta;
  delete[] xVars;
  delete[] yVars;
  d_mp.update();
  d_sub.update();
}



















