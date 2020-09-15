#include "master.h"


bool Master::add_cut(BendersCut cut, Solution sol, double tol)
{
      // adding the constraint: (1 + tau) theta + beta^T x >= alpha  
  double alpha_betax = -inner_product(cut.d_beta.begin(), cut.d_beta.end(), sol.xVals.begin(), -cut.d_alpha);    // alpha - beta^T x
  double kappa = 1 + cut.d_tau;
  //cout << "kappa = " << cut.d_tau << '\n';
  //cout << "old theta = " << sol.thetaVal << ". new theta = " << alpha_betax - cut.d_tau * sol.thetaVal << ".\n";
  //cout << "theta = " << alpha_betax / (1 + cut.d_tau) << '\n';

  bool add_cut = alpha_betax > kappa * sol.thetaVal + tol;
  
  if (add_cut) // then add cut and return false
  {
    ++d_nSlacks;

        // adding the cut to d_cmodel
    GRBaddvar(d_cmodel, 0, NULL, NULL, 0, 0, GRB_INFINITY, GRB_CONTINUOUS, NULL);  // slack
    
    size_t numVars = d_n1 + 2; // theta, x-vars, and slack  
    
    int cind[numVars];
    iota(cind, cind + d_n1 + 1, 0);
    cind[d_n1 + 1] = d_n1 + d_nSlacks;
    
    double cval[numVars];
    cval[0] = kappa;
    copy_n(cut.d_beta.begin(), d_n1, cval + 1);

    cval[numVars - 1] = -1;     // >= constraint, so slack features with -1
    GRBaddconstr(d_cmodel, numVars, cind, cval, GRB_EQUAL, cut.d_alpha - kappa * d_L, NULL);
    
    // update slack identities
    d_kappa.push_back(kappa);
    d_beta.push_back(cut.d_beta);
    d_gamma.push_back(cut.d_alpha);

    GRBupdatemodel(d_cmodel);

    GRBLinExpr lhs = kappa * d_theta;
    lhs.addTerms(cut.d_beta.data(), d_xvars.data(), d_xvars.size());
    d_interceptor.addConstr(lhs >= cut.d_alpha);
    d_interceptor.update();


    return false;
  }
  return true; // betaxgamma >= theta, no cut added 
}