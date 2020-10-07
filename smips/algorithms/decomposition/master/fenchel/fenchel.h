#ifndef TSSMIPS_FENCHEL_H
#define TSSMIPS_FENCHEL_H

#include "gurobi_c++.h"
#include <vector>

#include "../../../../problem_data/problem.h"
#include "../../cut/benderscut.h"

class Fenchel
{
  public:
    Problem &d_problem;
      // cuts are of the form kappa * theta + beta^T x >= alpha
    GRBModel d_mp;
    GRBVar d_alpha;
    vector<GRBVar> d_beta;
    GRBVar d_kappa;

    GRBModel d_sub;
    vector<GRBVar> d_xvars;
    GRBVar d_theta;

    struct Point
    {
        vector<double> d_xvals;
        double d_theta;
        double d_lb;
        double d_ub;
    };

    vector<Point> d_points;
    bool d_rcut;

    Fenchel(Problem &problem, GRBEnv &env);
    Fenchel(const Fenchel &other);

    void update_bounds(int var, double val, bool lower);
    void reverse_cut(double UB);
    void add_row(BendersCut &cut);

    BendersCut fenchel_cut(vector<double> &x, double theta, double tol, bool reset = false);
    bool solve_mp(double tol);
    BendersCut get_candidate();
    Point solve_sub();
    void set_mp_obj(vector<double> &x, double &theta);  // takes (x, theta) and sets master objective coefficients
    void set_sub_obj(BendersCut &cut);                        // takes (alpha, beta, tau) and sets subproblem coefficients
    void add_mp_cut(Point const &point);

    void clear_mp();
    void chg_mp_tol(bool focus);
};

#endif //TSSMIPS_FENCHEL_H
