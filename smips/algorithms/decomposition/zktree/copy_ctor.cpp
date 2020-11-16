#include "zktree.h"

ZkTree::ZkTree(const ZkTree &other)
:
  d_problem(other.d_problem),
  d_cglp(other.d_cglp),
  d_L(other.d_L),
  d_lb_mult_inds(other.d_lb_mult_inds),
  d_ub_mult_inds(other.d_ub_mult_inds),
  d_rcut_inds(other.d_rcut_inds),
  d_cglp_val(other.d_cglp_val),
  d_candidate(other.d_candidate)
{
  d_nodes.reserve(other.d_nodes.size());
  for (ZK *node : other.d_nodes)
  {
    ZK *node_copy = new ZK(*node);
    d_nodes.push_back(node_copy);
  }

  d_alpha = d_cglp.getVarByName("alpha");
  d_tau = d_cglp.getVarByName("tau");
  d_kappa = d_cglp.getVarByName("kappa");

  size_t n1 = other.d_beta.size();

  d_beta.reserve(n1);
  for (size_t var = 0; var != n1; ++var)
    d_beta.push_back(d_cglp.getVarByName("beta_" + to_string(var)));

  size_t nTerms = d_nodes.size();
  d_lambda.reserve(nTerms);
  for (size_t term = 0; term < nTerms; ++term)
  {
    vector<GRBVar> mults;
    size_t nMults = other.d_lambda[term].size();
    mults.reserve(nMults);
    for (size_t mult = 0; mult != nMults; ++mult)
      mults.push_back(d_cglp.getVarByName("lambda_" + to_string(term) + "_" + to_string(mult)));

    d_lambda.push_back(mults);
  }


  GRBConstr *con_ptr = d_cglp.getConstrs();
  size_t con_idx = 0;
  size_t nCons = n1 + 3;
  d_constrs.reserve(nTerms);
  for (size_t term = 0; term != nTerms; ++term)
  {
    vector<GRBConstr> constrs;
    constrs.reserve(nCons);
    for (size_t con = 0; con != nCons; ++con)
    {
      constrs.push_back(con_ptr[con_idx]);
      ++con_idx;
    }
    d_constrs.push_back(constrs);
  }
  delete[] con_ptr;

  d_cglp.update();
}

