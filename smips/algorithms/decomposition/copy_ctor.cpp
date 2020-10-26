#include "benders.h"

Benders::Benders(const Benders &other)
:
  d_problem(other.d_problem),
  d_env(other.d_env),
  d_n1(other.d_n1), 
  d_p1(other.d_p1),
  d_m2(other.d_m2), 
  d_n2(other.d_n2), 
  d_S(other.d_S),
  d_master(other.d_master),
  d_sub(other.d_sub),
  d_lr(other.d_lr),
  d_gomory(other.d_gomory),
  d_ald(other.d_ald),
  d_pslp(other.d_pslp),
  d_agg(other.d_agg),
  d_def(other.d_def),
  d_lb(other.d_lb),
  d_ub(other.d_ub),
  d_visited(other.d_visited),
  d_objectives(other.d_objectives),
  d_UB(other.d_UB)
{
  d_xvals = new double[d_n1];
  d_incumbent = new double[d_n1];
  
  copy_n(other.d_xvals, d_n1, d_xvals);
  copy_n(other.d_incumbent, d_n1, d_incumbent);
}