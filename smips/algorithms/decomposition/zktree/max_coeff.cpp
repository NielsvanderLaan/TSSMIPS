#include "zktree.h"

static bool abs_compare(double a, double b)
{
  return (std::abs(a) < std::abs(b));
}

double ZkTree::max_coeff(BendersCut &cut)
{
  double beta_max = abs(*max_element(cut.d_beta.begin(), cut.d_beta.end(), abs_compare));

  return max(beta_max, max(abs(cut.d_alpha), abs(cut.d_tau)));
}