
#include <mex.h>

#include <math.h>
#include <iostream>
#include <vector>
#include <set>
#include <eigen3/Eigen/Dense>
#include <list>
#include "fastQP.h"
#include <eigen3/Eigen/LU>



using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

double log_expr(double input)
{
    return log(input);
}

double invert_expr(double input)
{
    return 1/input;
}

int NewtonIP(MatrixXd& Q, MatrixXd& A, MatrixXd& b, VectorXd& c, int T, VectorXd& x0, int maxIter, double threshold,
        double alpha, double beta, int newtonIter, double f)
{
    int feasible;
    cout << "got here" << endl;
    cout << A << endl;
    cout << x0 << endl;
    int m = A.rows();
    int n = A.cols();
    //consider changing it to single vector
    MatrixXd x(m,n);
    x.col(0) = x0;
     VectorXd dx(m);
     MatrixXd hess(m,m);
//     //note that this doesn't have to be done every time.
     int q = -10;
     double s0 = -1.0/q;
     MatrixXd C(3,3);
     C << 2.0,0.0,0.0,
          2*s0,1.0,0.0,
          0.0,0.0,1.0;
     VectorXd r(3);
     r << 1.0/pow(s0,2.0),
         -1.0/s0,
         -log(s0);
     VectorXd a = C.lu().solve(r);
     
    cout << "how about here" << endl;
     
    //newtonIter
     for (int i = 0; i < newtonIter; i++)
     {
         VectorXd xi = x.col(i);
         VectorXd s = A*xi - b;
         cout << s << endl;
         cout << "strange as" << endl;
         //1 for log evaluation 0 for quadratic
         VectorXd d;
         VectorXd grad;
         MatrixXd Hess;
         MatrixXd Afp;
         MatrixXd Aup;
         VectorXd bfp;
         VectorXd bup;
         int q_eval;
         int l_eval;
         
         int quad = 0;
         for (int k = 0; k < m; k++)
         {
             if (s(k) <= 0)
             {
                quad = 1;
            }
         }
         
         
         if (quad)
         {
             cout << "shouldn't be here" << endl;
             vector<int> l;
             vector<int> q;
             for (int k = 0; k < m; k++)
             {
                 if (s(k) <= s0)
                 {
                     l.push_back(k);
                 } else
                 {
                     q.push_back(k);
                 }
             }
             q_eval = q.size();
             l_eval = l.size();
             
             MatrixXd Af(q_eval,n);
             MatrixXd Au(l_eval,n);
             VectorXd bf(q_eval);
             VectorXd bu(l_eval);
             
             for (int k = 0; k < q_eval; k++)
             {
                 Af.row(k) = A.row(l[k]);
                 bf(k) = b(l[k]);
             }
             for (int k = 0; k < l_eval; k++)
             {
                 Au.row(k) = A.row(q[k]);
                 bu(k) = b(q[k]);
             }
             VectorXd fexpr = Af*xi - bf;
             VectorXd uexpr = Au*xi - bu;
             MatrixXd expr;
             expr << uexpr.array().square(), uexpr, VectorXd::LinSpaced(l_eval,1,1);
 

             VectorXd uquad = expr*a;
             double fe = fexpr.array().log().sum();
             double un = uquad.sum();
             f = -fe + un;
             d = uexpr.array().inverse();
             grad = -Af.transpose()*d + (2*a(0)*uexpr.asDiagonal()*Au + a(1)*Au).transpose()*VectorXd::LinSpaced(l_eval,1,1);
             VectorXd sub = d.array().square();
             Hess = Af.transpose()*sub.asDiagonal()*Af + 2*a(0)*(Au.transpose()*Au); 
             Afp = Af;
             Aup = Au;
             bfp = bf;
             bup = bu;
         } else
         {
             //VectorXd slog = s.array().log();
             f = (double)(T*xi.transpose()*Q*xi) + (double)(T*c.transpose()*xi) - ((VectorXd) s.array().log()).sum();
             d = s.array().inverse();
             grad = T*(Q + Q.transpose())*xi + T*c - A.transpose()*d;
             Hess = T*(Q + Q.transpose()) + A.transpose() * (VectorXd)(d.array().square())*A;
         }
         //Newton Step and Decrement
         VectorXd Dx = - Hess.lu().solve(grad);
         double lambda = - grad.transpose() * Dx;
         
         //Stopping Criterion
         int exit = 0;
         for (int k = 0; k < n; k++)
         {
             if (std::isinf((double) Dx(k)) || std::isnan((double) Dx(k))) exit = 1;
         }
         if (lambda/abs(f) < threshold || i == newtonIter || exit)
         {
            if (lambda/abs(f) < threshold)
            {
                cout << "threshold exit" << endl;
            } else if (i == newtonIter)
            {
                cout << "iteration exit" << endl;
            } else
            {
                cout << "boundary collision exit" << endl;
            }
            if (quad)
            {
                feasible = 0;
            } else
            {
                feasible = 1;
            }
            break;
         } else 
         {
             //backtracking line search
             double t = 1.0;
             if (!quad)
             {
                VectorXd bound = b - A*xi;
                VectorXd adx = A*Dx;
                for (int k = 0; k < m; k++)
                {
                    double slack = bound(k)/adx(k);
                    if (!std::isnan(slack) && !std::isinf(slack))
                    {
                        if (slack > 0)
                        {
                            t = min(t,slack);
                        }
                    }
                }
                t = t*.999;
             }
             int backtrack = 1;
             VectorXd x_new;
             double f_new;
             while (backtrack)
             {
                 x_new = xi + t*Dx;
                 if (quad)
                 {
                     VectorXd fexpr = Afp*x_new - bfp;
                     VectorXd uexpr = Aup*x_new - bup;
                     MatrixXd expr;
                     expr << uexpr.array().square(), uexpr, VectorXd::LinSpaced(l_eval,1,1);
                     VectorXd uquad = expr*a;
                     double fe = fexpr.array().log().sum();
                     double un = uquad.sum();
                     f_new = -fe + un;

                 } else 
                 {
                     VectorXd s = A*x_new - b;
                     f_new = (double)(T*x_new.transpose()*Q*x_new) + (double)(T*c.transpose()*x_new) - ((VectorXd) s.array().log()).sum();  
                 }
                 if (f_new < f + alpha*t*grad.transpose()*Dx)
                 {
                     backtrack = 0;
                 } else
                 {
                     t = beta*t;
                 }
             }
             x.col(i+1) = x_new;
//              if (std::isnan(x_new) || std::isinf(x_new))
//              {
//                  cout << "holy crap" << endl;
//              }
         }
     }
    return feasible;
}

