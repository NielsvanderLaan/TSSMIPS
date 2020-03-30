#include "master.h"

Master::Master(GRBEnv &env, GRBenv *c_env, Problem &problem, bool zk_safe)
{
  size_t n1 = problem.d_n1;
  size_t m1 = problem.d_m1;
  size_t p1 = problem.d_p1;
  
  size_t fs_leq, fs_geq;
  fs_leq = problem.d_fs_leq; fs_geq = problem.d_fs_geq; 
 
  d_n1 = n1; 
  d_p1 = p1;
  d_L = problem.d_L;
  d_nSlacks = fs_leq + fs_geq;
  
  vector<vector<double>> &Amat = problem.d_Amat;
  
  double *c = problem.d_c.data();
  double *rhs = problem.d_b.data();  
  double *l1 = problem.d_l1.data();
  double *u1 = problem.d_u1.data();
  
  // instantiating c-api gurobi model (in order to use advanced simplex routines)
  
  GRBnewmodel(c_env, &d_cmodel, NULL, 0, NULL, NULL, NULL, NULL, NULL); 
  GRBaddvar(d_cmodel, 0, NULL, NULL, 1.0, problem.d_L, GRB_INFINITY, GRB_CONTINUOUS, NULL);     // theta
  
  if (zk_safe)
  {  
      // we do not impose integer requirements on xvars: we solve the IP manually via B&B
      // we do not impose upper and lower bounds directly: we include them as constraints (to enable the B&B and to compute the simplex tableau)
    GRBaddvars(d_cmodel, n1, 0, NULL, NULL, NULL, c, NULL, NULL, NULL, NULL);  // xvars 
    
    double val = 1;
    for (size_t idx = 0; idx != n1; ++idx)
    {
      int var = idx + 1;
      GRBaddconstr(d_cmodel, 1, &var, &val, GRB_EQUAL, l1[idx], NULL);
    }
    for (size_t idx = 0; idx != n1; ++idx)
    {
      int var = idx + 1;
      GRBaddconstr(d_cmodel, 1, &var, &val, GRB_EQUAL, u1[idx], NULL);
    }
    
    // adding slacks (x - s = l (n1 times), x + s = u (n1 times))
    size_t nBounds = 2 * n1;
    int vbeg[nBounds];
    iota(vbeg, vbeg + nBounds, 0);
    int *vind = vbeg;
    double vval[nBounds];
    fill_n(vval, n1, -1);
    fill_n(vval + n1, n1, 1);
    GRBaddvars(d_cmodel, nBounds, nBounds, vbeg, vind, vval, NULL, NULL, NULL, NULL, NULL);   
    
    // storing slack identities 
    for (size_t var = 0; var != n1; ++var)
    {
      d_kappa.push_back(0);
      vector<double> row(n1);
      row[var] = 1;
      d_beta.push_back(row);
      d_gamma.push_back(l1[var]);
    }
    for (size_t var = 0; var != n1; ++var)
    {
      d_kappa.push_back(0);
      vector<double> row(n1);
      row[var] = -1;
      d_beta.push_back(row);
      d_gamma.push_back(-u1[var]);
    }
    
    d_nSlacks += 2 * n1;
  } else
  {
    char vtypes[n1];
    fill_n(vtypes, p1, GRB_INTEGER);
    fill(vtypes + p1, vtypes + n1, GRB_CONTINUOUS);
    GRBaddvars(d_cmodel, n1, 0, NULL, NULL, NULL, c, l1, u1, vtypes, NULL);  // xvars
  } 
  
      // adding constraints
  int cind[n1];
  iota(cind, cind + n1, 1);
  
  for (size_t con = 0; con != m1; ++con)
    GRBaddconstr(d_cmodel, n1, cind, Amat[con].data(), GRB_EQUAL, rhs[con], NULL);
    
  
      // adding slacks
  size_t con_idx = zk_safe ? 2 * n1 : 0;      // constraint index corresponding to first fs constraint  
           
  size_t nSlacks = fs_leq + fs_geq;
  int vbeg[nSlacks];
  iota(vbeg, vbeg + nSlacks, 0);
  int vind[nSlacks];
  iota(vind, vind + nSlacks, con_idx);        // con_idx used here
  double vval[nSlacks];
  fill_n(vval, fs_leq, 1);
  fill_n(vval + fs_leq, fs_geq, -1);
  GRBaddvars(d_cmodel, nSlacks, nSlacks, vbeg, vind, vval, NULL, NULL, NULL, NULL, NULL);   
  
      // storing slack identities  
  for (size_t con = 0; con != fs_leq; ++con)
  {
    d_kappa.push_back(0);
    vector<double> minus_beta(Amat[con]);
    for_each(minus_beta.begin(), minus_beta.end(), [](double &val){ val *= -1.0; });
  
    d_beta.push_back(minus_beta);
    d_gamma.push_back(-rhs[con]); 
  }  
  
  for (size_t con = fs_leq; con != fs_leq + fs_geq; ++con)
  {
    d_kappa.push_back(0);
    d_beta.push_back(Amat[con]);
    d_gamma.push_back(rhs[con]); 
  }  
  
  GRBupdatemodel(d_cmodel);
  //GRBwrite(d_cmodel, "master.lp");
}








