#include "zk.h"

ZK::ZK(GRBenv *env, GRBEnv &cpp_env, Problem &problem, size_t scenario, bool lap)
:
  d_n1(problem.d_n1), 
  d_p1(problem.d_p1), 
  d_n2(problem.d_n2), 
  d_p2(problem.d_p2), 
  d_m2(problem.d_m2), 
  d_nConstrs(d_m2),
  d_cglp(problem, cpp_env, scenario, lap),
  d_Tmat(problem.d_Tmat),
  d_tau(d_m2),
  d_omega(problem.d_omega[scenario]),
  d_signs(d_m2),
  d_lb_inds(d_n2, -1),
  d_ub_inds(d_n2, -1),
  d_L(problem.d_L),
  d_yvals(d_n2)
{
  d_Wmat = problem.d_fix_rec ? problem.d_Wmat : problem.d_W_omega[scenario];

  GRBnewmodel(env, &d_model, NULL, 0, NULL, NULL, NULL, NULL, NULL);  
  GRBsetintparam(GRBgetenv(d_model), "ScaleFlag", 0);
  GRBsetintparam(GRBgetenv(d_model), "Method", 1);

  // adding variables  (assumed continuous, lower and upper bounds imposed later in canonical form)
  vector<double> &q = problem.d_fix_rec ? problem.d_q : problem.d_q_omega[scenario];
  GRBaddvars(d_model, d_n2, 0, NULL, NULL, NULL, q.data(), NULL, NULL, NULL, NULL);
    
  size_t nLeq = problem.d_ss_leq;
  size_t nGeq = problem.d_ss_geq;
        // updating signs of inequalities (default = 0)
  fill_n(d_signs.begin(), nLeq, -1);        // <= constraints
  fill_n(d_signs.begin() + nLeq, nGeq, 1);  // >= contraints 

  int nSlacks = nLeq + nGeq;        // number of slack variables (c.t. <= and >= constraints)

  vector<double> vval(nSlacks);
  fill_n(vval.begin(), nLeq, 1.0);          // coefficient of slack variable (c.t. <= constraints)
  fill_n(vval.begin() + nLeq, nGeq, -1.0);  // idem (c.t. >= constraints)
  
  for (size_t var = 0; var != d_n2; ++var)    // looping over variables: checking for non-trivial upper and lower bounds
  {
    double lb = problem.d_l2[var];
    double ub = problem.d_u2[var];
    if (lb > 0)                                // non-trivial lb
    {
      d_Wmat.push_back(vector<double>(d_n2));  // update constraint data
      d_Wmat[d_Wmat.size() - 1][var] = 1;      // idem
      d_Tmat.push_back(vector<double>(d_n1));  // idem
      d_tau.push_back(0);                      // idem
      d_omega.push_back(lb);                   // idem 
      d_signs.push_back(1);                    // >= constraint
      vval.push_back(-1.0);                    // coefficient of slack variable  
      d_lb_inds[var] = d_nConstrs;             // index of constraint corresponding to bound (for updating purposes) 
      ++d_nConstrs;                            // we added the constraint x >= l
      ++nSlacks;                               // using one slack variable 
    }
    if (ub < 1e10)                             // idem 
    {
      d_Wmat.push_back(vector<double>(d_n2));
      d_Wmat[d_Wmat.size() - 1][var] = 1;
      d_Tmat.push_back(vector<double>(d_n1));
      d_tau.push_back(0);
      d_omega.push_back(ub);
      d_signs.push_back(-1);
      vval.push_back(1.0);
      d_ub_inds[var] = d_nConstrs;
      ++d_nConstrs;
      ++nSlacks;
    }
  }
  d_nVars = d_n2 + nSlacks;       // total number of variables 
  
        // adding constraints
  int cind[d_n2];                  // variable indices
  iota(cind, cind + d_n2, 0);      
  for (size_t con = 0; con != d_nConstrs; ++con)    
    GRBaddconstr(d_model, d_n2, cind, d_Wmat[con].data(), GRB_EQUAL, 0, NULL);    // no slacks (yet)
  
        // adding slacks  
  int vbeg[nSlacks];
  iota(vbeg, vbeg + nSlacks, 0);
  int vind[nSlacks];                               // constraint indices 
  iota(vind, vind + nLeq + nGeq, 0);               // corresponding to <= and >= constraints
  iota(vind + nLeq + nGeq, vind + nSlacks, d_m2);  // corresponding to bounds (x >= l and x <= u)
  GRBaddvars(d_model, nSlacks, nSlacks, vbeg, vind, vval.data(), NULL, NULL, NULL, NULL, NULL);

  GRBupdatemodel(d_model);
}
















