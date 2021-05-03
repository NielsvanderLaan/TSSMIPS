#include "problem.h"

Problem::Problem(Data &generator, GRBEnv &env)
:
d_fix_rec(true),
d_fix_tech(true),
d_gen(generator),
d_env(env),
d_L(0)
{}