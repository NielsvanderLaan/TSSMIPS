#include "problem.h"

Problem::Problem(Data &generator, GRBEnv &env)
:
d_gen(generator),
d_env(env),
d_L(0)
{}