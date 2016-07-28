
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


int NewtonIP(MatrixXd& Q, MatrixXd& A, MatrixXd& b, VectorXd& c, int T, VectorXd& x0, int maxIter, double threshold,
        double alpha, double beta, int newtonIter, double f)
{
    int feasible;
    //cout << "got here" << endl;
    //cout << A << endl;
    //cout << x0 << endl;
    int m = A.rows();
    int n = A.cols();
    //consider changing it to single vector
    MatrixXd x(n,newtonIter + 1);
    x.col(0) = x0;
    //VectorXd dx(m);
    //MatrixXd Hess(m,m);
//     //note that this doesn't have to be done every time.
    int q = -10;
    double s0 = -1.0/q;
    MatrixXd E(3,3);
    E << 2.0,0.0,0.0,
            2*s0,1.0,0.0,
            0.0,0.0,1.0;
    VectorXd r(3);
    r << 1.0/pow(s0,2.0),
            -1.0/s0,
            -log(s0);
    VectorXd a = E.colPivHouseholderQr().solve(r);
    
    int final_index;
    
    //cout << "how about here" << endl;
     
    //newtonIter
     for (int i = 0; i < newtonIter; i++)
     {
         feasible = 1;
         VectorXd xi = x.col(i);
         
         //cout << xi << endl;
         
         VectorXd s = A*xi - b;
    
 //         cout << "xi" << endl;
  //        cout << xi << endl;
//          cout << "s" << endl;
//          cout << s << endl;
      //   cout << "strange as" << endl;
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
      //   cout << "Hope for Best" << endl;
         
         if (quad)
         {
             feasible = 0;
             vector<int> l;
             vector<int> q;
             for (int k = 0; k < m; k++)
             {
                 if (s(k) <= s0)
                 {
                     q.push_back(k);
                 } else
                 {
                     l.push_back(k);
                 }
             }
             //cout << "DEATH" << endl;
             q_eval = q.size();
             l_eval = l.size();
             
             MatrixXd Af(l_eval,n);
             MatrixXd Au(q_eval,n);
             MatrixXd bf(l_eval,1);
             MatrixXd bu(q_eval,1);
             //cout << "TOO MUCH DEATH" << endl;
             //cout << q << endl;
             //cout << l << endl;
             
             for (int k = 0; k < l_eval; k++)
             {
                 Af.row(k) = A.row(l[k]);
                 bf.row(k) = b.row(l[k]);
             }
             //cout << "make it here" << endl;
             for (int k = 0; k < q_eval; k++)
             {
                 Au.row(k) = A.row(q[k]);
                 bu.row(k) = b.row(q[k]);
             }
             //cout << "CONFIRMED" << endl;
             VectorXd fexpr = Af*xi - bf;
             VectorXd uexpr = Au*xi - bu;
             MatrixXd expr(q_eval,3);
//              cout << "Probably die here" << endl;
//              cout << uexpr.array().square() << endl;
//              cout << uexpr << endl;
//              cout << VectorXd::LinSpaced(q_eval,1,1) << endl;
             expr << uexpr.array().square(), uexpr, VectorXd::LinSpaced(q_eval,1,1);
             //expr << uexpr.array().square(), uexpr;

             VectorXd uquad = expr*a;
             double fe = fexpr.array().log().sum();
             double un = uquad.sum();
             f = -fe + un;
             d = fexpr.array().inverse();
             grad = -Af.transpose()*d + (2*a(0)*uexpr.asDiagonal()*Au + a(1)*Au).transpose()*VectorXd::LinSpaced(q_eval,1,1);
             VectorXd sub = d.array().square();
             Hess = Af.transpose()*sub.asDiagonal()*Af + 2*a(0)*(Au.transpose()*Au); 
             Afp = Af;
             Aup = Au;
             bfp = bf;
             bup = bu;
         } else
         {
             //VectorXd slog = s.array().log();
      //       cout << xi << endl;
             f = (double)(T*xi.transpose()*Q*xi) + (double)(T*c.transpose()*xi) - ((VectorXd) s.array().log()).sum();
             d = s.array().inverse();
             //cout << d << endl;
             grad = T*(Q + Q.transpose())*xi + T*c - A.transpose()*d;
             Hess = T*(Q + Q.transpose()) + A.transpose() * ((VectorXd)(d.array().square())).asDiagonal()*A;
             //cout << Q << endl;
             //cout << A << endl;
       //      cout << "NO Problem Here" << endl;
         }
         //Newton Step and Decrement
         VectorXd Dx = - Hess.ldlt().solve(grad);
         double lambda = - grad.transpose() * Dx;
         
//           cout << "xi" << endl;
//           cout << xi << endl;
//           cout << "grad" << endl;
//           cout << grad << endl;
//           cout << "Hess" << endl;
//           cout << Hess << endl;
//           cout << "Dx" << endl;
//           cout << Dx << endl;
         
      //   cout << "stopping criterion" << endl;
         //Stopping Criterion
         int exit = 0;
         for (int k = 0; k < n; k++)
         {
             if (std::isinf((double) Dx(k)) || std::isnan((double) Dx(k))) exit = 1;
             // cout << "Dx FOR REAL" << endl;
//              cout << std::isinf((double) Dx(k)) << endl;
//              cout << std::isnan((double) Dx(k)) << endl;
             // cout << Dx << endl;
         }
         if (lambda/abs(f) < threshold || exit)
         {
            if (lambda/abs(f) < threshold)
            {
                //cout << "threshold exit" << endl;
            } else 
            {
                //cout << "boundary collision exit" << endl;
            }
            //bit of tricky off by ones
            final_index = i;
            break;
         } else 
         {
             //backtracking line search
             double t = 1.0;
             VectorXd bound; //distance from bound
             VectorXd adx; // step direction
             int fcon; //number of feasible constraints
             if (quad)
             {
                 bound = bfp - Afp*xi;
                 adx = Afp*Dx;
                 fcon = l_eval;
             } else
             {
                 bound = b - A*xi;
                 adx = A*Dx;
                 fcon = m;
             }
             for (int k = 0; k < fcon; k++)
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
             
//              if (!quad)
//              {
//                 VectorXd bound = b - A*xi;
//                 VectorXd adx = A*Dx;
//                 for (int k = 0; k < m; k++)
//                 {
//                     double slack = bound(k)/adx(k);
//                     if (!std::isnan(slack) && !std::isinf(slack))
//                     {
//                         if (slack > 0)
//                         {
//                             t = min(t,slack);
//                         }
//                     }
//                 }
//                 t = t*.999;
//              }
             int backtrack = 1;
             VectorXd x_new;
             double f_new;
             while (backtrack)
             {
                  //cout << "stuck" << endl;
                  x_new = xi + t*Dx;
//                  cout << "x_new" << endl;
//                  cout << x_new << endl;
//                  cout << "Dx" << endl;
//                  cout << Dx << endl;
                 if (quad)
                 {
                     //cout << "NOOOOOOOOOO" << endl;
                     VectorXd fexpr = Afp*x_new - bfp;
                     VectorXd uexpr = Aup*x_new - bup;
                     MatrixXd expr(q_eval,3);
                     expr << uexpr.array().square(), uexpr, VectorXd::LinSpaced(q_eval,1,1);
                     VectorXd uquad = expr*a;
                     //cout << "complex log" << endl;
                     //cout << complex_log(-2) << endl;
                     double fe = fexpr.array().log().sum();
                     double un = uquad.sum();
                     //cout << "fe" << endl;
                     //cout << fe << endl;
                     //cout << "un" << endl;
                     //cout << un << endl;
                     f_new = -fe + un;

                 } else 
                 {
                     VectorXd s = A*x_new - b;
                     f_new = (double)(T*x_new.transpose()*Q*x_new) + (double)(T*c.transpose()*x_new) - ((VectorXd) s.array().log()).sum();  
                 }
//                  cout << "s" << endl;
//                  cout << s << endl;
//                  cout << "log terms" << endl;
//                  cout << ((VectorXd) s.array().log()).sum() << endl;
//                   cout << "Q" << endl;
//                   cout << Q << endl;
//                   cout << "c" << endl;
//                   cout << c << endl;
//                   cout << T << endl;
//                   cout << (double)(T*x_new.transpose()*Q*x_new) << endl;
//                   cout << (double)(T*c.transpose()*x_new) << endl;
//                 
//                  cout << "innocent terms"  << endl;
//                  cout << (double)(T*x_new.transpose()*Q*x_new) + (double)(T*c.transpose()*x_new) << endl;
                  //cout << "x_new" << endl;
                  //cout << x_new << endl;
                  //cout << "f_new" << endl;
                  //cout << f_new << endl;
//                  cout << "f" << endl;
//                  cout << f << endl;
//                  cout << "linear f" << endl;
//                  cout << f + alpha*t*grad.transpose()*Dx << endl;
                 
                 if (f_new < f + alpha*t*grad.transpose()*Dx)
                 {
                     backtrack = 0;
                 } else
                 {
                     t = beta*t;
                 }
                 if (t < .00001)
                 {
                     break;
                 }
             }
             //cout << "unstuck" << endl;
             x.col(i+1) = x_new;
             final_index = i+1;
         }
     }
    x0 = x.col(final_index);
    return feasible;
}

