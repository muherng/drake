#include <mex.h>

#include <math.h>
#include <iostream>
#include <vector>
#include <set>
#include <eigen3/Eigen/Dense>
#include "fastQP.h"
#include "fastIP.cpp"

using namespace Eigen;
using namespace std;


/*
 * timestwo.c - example found in API guide
 *
 * Computational function that takes a scalar and doubles it.
 *
 * This is a MEX-file for MATLAB.
 * Copyright 1984-2011 The MathWorks, Inc.
 */


void timestwo(double y[], double x[])
{
    y[0] = 2.0*x[0];
}

void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[] )
{
    //my own code
    if (nrhs < 1) {
        mexErrMsgIdAndTxt(
                "Drake:fastIP:NotEnoughInputs",
                "Usage [x, info, active] = fastQP(Q, f[, Aeq, beq, Ain, bin, active])");
    }
    int arg = 0, nblks = 1;
    if (mxIsCell(prhs[arg])) nblks = mxGetNumberOfElements(prhs[arg]);
    MatrixXd* Q = new MatrixXd[nblks];
    vector<MatrixXd*> QblkMat;
    
    if (mxIsCell(prhs[arg])) {
        mxArray* QblkDiagCellArray = (mxArray*)prhs[arg++];
        for (int i = 0; i < nblks; i++) {
            mxArray* Qblk = mxGetCell(QblkDiagCellArray, i);
            int m = mxGetM(Qblk), n = mxGetN(Qblk);
            if (m * n == 0)
                continue;
            else if (m == 1 || n == 1)  // then it's a vector
                Q[i] = Map<MatrixXd>(mxGetPr(Qblk), m * n, 1);
            else
                Q[i] = Map<MatrixXd>(mxGetPr(Qblk), m, n);
            QblkMat.push_back(&Q[i]);
        }
        if (QblkMat.size() < 1)
            mexErrMsgIdAndTxt("Drake:FastQP:BadInputs", "Q is empty");
    } else {
        int m = mxGetM(prhs[arg]), n = mxGetN(prhs[arg]);
        if (m * n == 0)
            mexErrMsgIdAndTxt("Drake:FastQP:BadInputs", "Q is empty");
        else if (m == 1 || n == 1)  // then it's a vector
            Q[0] = Map<MatrixXd>(mxGetPr(prhs[arg]), m * n,
                    1);  // always want a column vector
        else
            Q[0] = Map<MatrixXd>(mxGetPr(prhs[arg]), m, n);
        arg++;
        QblkMat.push_back(&Q[0]);
    }
    
    int m = 0;
    int n = 0;
    
    VectorXd* c = new VectorXd[mxGetNumberOfElements(prhs[arg])];
    if (nrhs > arg && mxGetNumberOfElements(prhs[arg]) > 0)
    {
        m = mxGetM(prhs[arg]);
        n = mxGetN(prhs[arg]);
        c[0] = Map<MatrixXd>(mxGetPr(prhs[arg]), m * n, 1);
        cout << c[0] << endl;
        arg++;
    }
     
    MatrixXd* A = new MatrixXd[mxGetNumberOfElements(prhs[arg])];
    if (nrhs > arg && mxGetNumberOfElements(prhs[arg]) > 0)
    {
        m = mxGetM(prhs[arg]);
        n = mxGetN(prhs[arg]);
        A[0] = Map<MatrixXd>(mxGetPr(prhs[arg]), m, n);
        cout << A[0] << endl;
        arg++;
    }
    
    
    VectorXd* b = new VectorXd[mxGetNumberOfElements(prhs[arg])];
    if (nrhs > arg && mxGetNumberOfElements(prhs[arg]) > 0)
    {
        m = mxGetM(prhs[arg]);
        n = mxGetN(prhs[arg]);
        b[0] = Map<MatrixXd>(mxGetPr(prhs[arg]), m * n, 1);
        arg++;
    }
    
    VectorXd* ub = new VectorXd[mxGetNumberOfElements(prhs[arg])];
    if (nrhs > arg && mxGetNumberOfElements(prhs[arg]) > 0)
    {
        m = mxGetM(prhs[arg]);
        n = mxGetN(prhs[arg]);
        ub[0] = Map<MatrixXd>(mxGetPr(prhs[arg]), m * n, 1);
        arg++;
    }
    
    VectorXd* lb = new VectorXd[mxGetNumberOfElements(prhs[arg])];
    if (nrhs > arg && mxGetNumberOfElements(prhs[arg]) > 0)
    {
        m = mxGetM(prhs[arg]);
        n = mxGetN(prhs[arg]);
        lb[0] = Map<MatrixXd>(mxGetPr(prhs[arg]), m * n, 1);
        arg++;
    }
    
    
    VectorXd* x = new VectorXd[mxGetNumberOfElements(prhs[arg])];
    if (nrhs > arg && mxGetNumberOfElements(prhs[arg]) > 0)
    {  
        m = mxGetM(prhs[arg]);
        n = mxGetN(prhs[arg]);
        x[0] = Map<MatrixXd>(mxGetPr(prhs[arg]), m * n, 1);
        arg++;
    }
    
    
    plhs[0] = mxCreateDoubleMatrix(m,n,mxREAL);
    double* outMatrix = mxGetPr(plhs[0]);
    fastIP(Q[0],c[0],A[0],b[0],ub[0],lb[0],x[0],outMatrix);
    
}