#include "pslp.h"

Pslp::Pslp(const Pslp &other)
: 
d_n1(other.d_n1),
d_S(other.d_S),
d_probs(other.d_probs)
{
  d_zk.reserve(d_S);
  
  for (size_t s = 0; s != d_S; ++s)
    d_zk.push_back(other.d_zk[s]);
}