#include "master.h"
Master::Master(const Master &other)
:
d_problem(other.d_problem),
d_n1(other.d_n1),
d_p1(other.d_p1),
d_L(other.d_L),
d_nSlacks(other.d_nSlacks),
d_zk_safe(other.d_zk_safe),
d_rcut_idx(other.d_rcut_idx),
d_kappa(other.d_kappa),
d_beta(other.d_beta),
d_gamma(other.d_gamma),
d_lb_con_inds(other.d_lb_con_inds),
d_lb_slack_inds(other.d_lb_slack_inds),
d_ub_con_inds(other.d_ub_con_inds),
d_ub_slack_inds(other.d_ub_slack_inds),
d_points(other.d_points)
{
  GRBupdatemodel(other.d_cmodel);
  d_cmodel = GRBcopymodel(other.d_cmodel);
}