#include "zk.h"

ZK::ZK(const ZK &other)
:
  d_n1(other.d_n1), 
  d_p1(other.d_p1),
  d_n2(other.d_n2), 
  d_p2(other.d_p2),
  d_m2(other.d_m2),
  d_nVars(other.d_nVars),
  d_nConstrs(other.d_nConstrs),
  d_cglp(other.d_cglp),
  d_Wmat(other.d_Wmat),
  d_Tmat(other.d_Tmat),
  d_tau(other.d_tau),
  d_omega(other.d_omega),
  d_signs(other.d_signs),
  d_lb_inds(other.d_lb_inds),
  d_ub_inds(other.d_ub_inds),
  d_yvals(other.d_yvals)
{
  GRBupdatemodel(other.d_model);
  d_model = GRBcopymodel(other.d_model);
}