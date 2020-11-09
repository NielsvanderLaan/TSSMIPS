#include "master.h"

void Master::update_bounds(int var, double val, bool lower)
{
  d_fenchel.update_bounds(var, val, lower);
  if (d_zk_safe)
  {
    vector<int> &con_inds = lower ? d_lb_con_inds : d_ub_con_inds;
    vector<int> &slack_inds = lower ? d_lb_slack_inds : d_ub_slack_inds;
    int con_idx = con_inds[var];
    if (con_idx == -1)
    {
      slack_inds[var] = d_nSlacks;
      GRBupdatemodel(d_cmodel);
      GRBgetintattr(d_cmodel, "NumConstrs", &con_inds[var]);  // set con idx of bound constraint to correct value
      GRBaddvar(d_cmodel, 0, NULL, NULL, 0.0, 0.0, GRB_INFINITY, GRB_CONTINUOUS, NULL);

      int inds[2] = {var + 1, 1 + d_n1 + d_nSlacks};     // var + 1 refers to the xvar (theta is skipped). 1 + d_n1 + d_nSlacks refers to slack variable
      double vals[2] = {1.0, lower ? -1.0 : 1.0};
      GRBaddconstr(d_cmodel, 2, inds, vals, GRB_EQUAL, val, NULL);    // x - slack = l1
      ++d_nSlacks;
          // updating slack variable identities
      d_kappa.push_back(0.0);
      vector<double> row(d_n1);
      row[var] = lower ? 1.0 : -1.0;
      d_beta.push_back(row);
      d_gamma.push_back(lower ? val : -val);
    } else
    {
      GRBsetdblattrelement(d_cmodel, "RHS", con_idx, val);
      d_gamma[slack_inds[var]] = lower ? val : -val;
    }
  } else  // not zk-safe
  {
    ++var;
    if (lower)
      GRBsetdblattrelement(d_cmodel, "LB", var, val);
    else
      GRBsetdblattrelement(d_cmodel, "UB", var, val);
  }
  GRBupdatemodel(d_cmodel);
}