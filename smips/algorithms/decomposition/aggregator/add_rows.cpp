#include "aggregator.h"

void Aggregator::add_rows(BendersCut &cut)
{
  for (size_t s = 0; s != d_cgmips.size(); ++s)
    d_cgmips[s].add_row(cut);
}