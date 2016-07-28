
#include <mex.h>

#include <math.h>
#include <iostream>
#include <vector>
#include <set>
#include <eigen3/Eigen/Dense>
#include "fastQP.h"
#include "NewtonIP.cpp"

using Eigen::MatrixXd;
using Eigen::VectorXd;

using namespace std;


void fastIP(MatrixXd& Q, VectorXd& c, MatrixXd& A, VectorXd& b, VectorXd& ub, VectorXd& lb, VectorXd& x, double* output)
{   
//Newton's method parameters
    int num_params = A.cols();
    MatrixXd A_full(A.rows() + 2*num_params,A.cols());
    A_full << A,
            MatrixXd::Identity(num_params,num_params),
            -1*MatrixXd::Identity(num_params,num_params);
    MatrixXd b_full(b.rows() + lb.rows() + ub.rows(), b.cols());
    b_full << b,
            lb,
            -1*ub;
   // cout << "display A full" << endl;
   // cout << A_full << endl;
    int maxIter = 3; // was 3
    int newtonIter = 3; //was 5
    double alpha = 0.01;
    double beta = 0.5;
            
     // Interior points parameters
     int T = 1; //inverse of log barrier constant
     double thresholdIP = 1*10^-5; //threshold for interior point method
     int mu = 100; //used to be 100  %schedule of T (consider revising this)
     double re = 1.0;
     //dualityGap = zeros(1, maxIter) ; %track duality gap
     //T_values = zeros(1, maxIter) ; %track values of T parameter
     int feasible = 0;
     double f = 0;
     
// Solve the program using Interior Point method with Newton + backtracking at each step
     int m = x.rows();
     int disp = 0;
     for (int j = 0; j < maxIter; j++)
     {
         feasible = NewtonIP(Q,A_full,b_full,c,T,x,maxIter,thresholdIP,alpha,beta,newtonIter,f);
         if (feasible)
         {
             T = mu*T;
         } else
         {
             T = (1/re)*T;
             disp = 1;
         }
         //cout << "SOLUTION TO ITERATION" << endl;
         if (disp)
         {
            //cout << x << endl;
         }
         //x0 = x.block(0,j,m,1);
     }
     
     for (int i = 0; i < x.rows(); i++)
     {
         output[i] = x(i);
     }
}