#ifndef CGMIP_H
#define CGMIP_H

#include <numeric>
#include "gurobi_c++.h"
#include "../../cut/benderscut.h"
#include "../../../../problem_data/problem.h"
#include "../../../../../debug.h"

using namespace std;

class CGMip
{
    // generate cuts of the form 
    // eta + beta^T x + tau theta >= alpha
  public: 
        // master problem
    Problem &d_problem;
    GRBModel d_mp;          
    GRBVar d_alpha;
    vector<GRBVar> d_beta;
    GRBVar d_tau;
        // sub problem
    GRBModel d_sub;         // d_sub is used to generate cut coefficients for d_mp
    vector<GRBVar> d_xVars;
    GRBVar d_theta;
    GRBVar d_eta;
    vector<GRBVar> d_yVars;
         
    struct Point
    {
      vector<double> d_x;
      double d_theta;
      double d_eta;
      double d_rhs_lb;
      double d_rhs_ub;
    };
    
    vector<Point> d_points;
    bool d_rcut;
 
    CGMip(GRBEnv &env, Problem &problem, size_t s);
    CGMip(const CGMip &other);
    CGMip(CGMip &&other) = delete;

    BendersCut generate_cut(double *x, double theta, bool init, double vwx, bool affine, double tol, bool int_feas, double &gap, double &npoints, bool reset = false);  // uses benders decomposition to find best cut
            // auxiliary functions for generate_cut()
    bool solve_mp(bool focus, bool affine, double M = 1e8);
    void update_mp();
    BendersCut get_candidate();
    Point solve_sub(bool focus = false);
    vector<Point> get_points();
    void set_mp_obj(double *x, double &theta);  // takes (x, theta) and sets master objective coefficients
    void set_sub_obj(BendersCut &cut);                        // takes (alpha, beta, tau) and sets subproblem coefficients
    void add_mp_cut(Point const &point);
    void clear_mp();
    
    void add_row(BendersCut &cut); // add the Benders' cut kappa theta >= beta^T x + gamma to d_sub
    void reverse_cut(double UB);
    void update_bound(size_t var, double val, bool lower);
    double mp_val();
    bool mp_optimal();
    bool check_mp_violation(double tol);

    double mp_max_coeff();
    void set_mp_bounds(double M, bool affine);

    double distance(Point &first, Point &second);
};


#endif





