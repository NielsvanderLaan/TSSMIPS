#include "zk.h"

void ZK::clear()
{
  GRBdelconstrs(d_model, d_cp_inds.size(), d_cp_inds.data());
  GRBdelvars(d_model, d_cp_slack_inds.size(), d_cp_slack_inds.data());
  GRBupdatemodel(d_model);

  d_nConstrs -= d_cp_inds.size();
  d_nVars -= d_cp_inds.size();

  for (size_t con = d_cp_inds.size() - 1; con != -1; --con)
  {
    size_t idx = d_cp_inds[con];
    d_Wmat.erase(d_Wmat.begin() + idx);
    d_Tmat.erase(d_Tmat.begin() + idx);
    d_tau.erase(d_tau.begin() + idx);
    d_omega.erase(d_omega.begin() + idx);
    d_signs.erase(d_signs.begin() + idx);
  }

  for (size_t var = 0; var != d_lb_inds.size(); ++var)
  {
    int idx = d_lb_inds[var];
    d_lb_inds[var] -= count_if(d_cp_inds.begin(), d_cp_inds.end(), [idx](int val){return val < idx;});

    idx = d_ub_inds[var];
    d_ub_inds[var] -= count_if(d_cp_inds.begin(), d_cp_inds.end(), [idx](int val){return val < idx;});
  }

  d_cp_inds.clear();
  d_cp_slack_inds.clear();
}