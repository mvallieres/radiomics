/*********************************************************************
 *     __          _           _    _  _____ 
 *   / _|        | |     /\  | |  | |/ ____|
 *  | |_ __ _ ___| |_   /  \ | |  | | |     
 *  |  _/ _` / __| __| / /\ \| |  | | |     
 *  | || (_| \__ \ |_ / ____ \ |__| | |____ 
 *  |_| \__,_|___/\__/_/    \_\____/ \_____|                                        
 *
 *  Usage:
 *   scores   = [10,5,90,1,-20,-1];
 *   labels   = [1,1,1,1,1,-1];
 *   posclass = 1;
 *   fastAUC(labels,scores,posclass)
 *
 *  fastAUC.cpp
 *      author  :   Enric Junqu� de Fortuny
 *      contact :   Applied Data Mining
 *                  Prinsstraat 13
 *                  University Antwerp
 *                  2000 Antwerp - Belgium
 *
 *                  <enric.junquedefortuny@ua.ac.be> 
 *                  http://ciri.be/
 *
 *      based on : "An introduction to ROC analysis" - Tom Fawcett, 
 *                 Pattern Recognition Letters 2005 
 * 21/06/2013   Bugfix abs/fabs (thank you Moritz von Heimendahl)
 ********************************************************************/
#include <matrix.h>
#include <mex.h>
#include <cmath>  // Added by Martin Vallières on June 7, 2015 (to use std::abs on floats)

#include <vector>
#include <algorithm>
#include <iostream>

/* Default mex definitions for compatibility and so on */
#ifndef MWSIZE_MAX
typedef int mwSize;
typedef int mwIndex;
typedef int mwSignedIndex;
#if (defined(_LP64) || defined(_WIN64)) && !defined(MX_COMPAT_32)
# define MWSIZE_MAX    281474976710655UL
# define MWINDEX_MAX   281474976710655UL
# define MWSINDEX_MAX  281474976710655L
# define MWSINDEX_MIN -281474976710655L
#else
# define MWSIZE_MAX    2147483647UL
# define MWINDEX_MAX   2147483647UL
# define MWSINDEX_MAX  2147483647L
# define MWSINDEX_MIN -2147483647L
#endif
#define MWSIZE_MIN    0UL
#define MWINDEX_MIN   0UL
#endif

double calcAUC      (double* labels, double* scores,int n,int posclass);
double trapezoidArea(double X1, double X2, double Y1, double Y2);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    /* Process arguments */
    if(nrhs < 2)
        mexErrMsgTxt("Usage: fastAUC(labels, scores [,posclass])");
    int poslabel = 1;
    if(nrhs > 2) {
        double * label  = mxGetPr(prhs[2]);
        poslabel        = label[0];
    }
    
    /* Associate inputs and output */
    double * labels = mxGetPr(prhs[0]);
    double * scores = mxGetPr(prhs[1]);    
    mxArray* ml_auc = plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
    double * auc    = mxGetPr(ml_auc );   

    /* Figure out dimensions and perform some checking */
    const mwSize* dims    = mxGetDimensions(prhs[0]);
    const mwSize* dims2   = mxGetDimensions(prhs[1]);
    int numdims     = mxGetNumberOfDimensions(prhs[0]);
    int dimy        = (int)dims[0]; 
    int dimx        = (int)dims[1];
    int numdims2    = mxGetNumberOfDimensions(prhs[1]);
    int dimy2       = (int)dims2[0]; 
    int dimx2       = (int)dims2[1];
    
    if( numdims != numdims2) {
            mexErrMsgTxt("Error: different dimensions for both given vectors.\n");
    } 
    else if( dimy != dimy2 || dimx != dimx2) {
           mexErrMsgTxt("Label and score vector should have the same size.\n");
    }
    
    /* The good stuff */
    auc[0] = calcAUC(labels,scores,std::max(dimx,dimy),poslabel);  
    return;
}

/* 
 * Initialize a pair list containing score/label pairs that need to be sorted to calculate the AUC 
 * Vector automatically sorts pair.first first and that is the desired behaviour.    
 *      scores      array of scores for each instance
 *      labels      array of labels for each instance
 *      n           amount of input instances
 *      posclass    label for the positive class
 *
 *      returns     the area under the ROC curve
 */
double calcAUC(double* labels, double * scores,int n,int posclass) {
	typedef std::pair<float,int> mypair;
	std::vector<mypair> L(n);
	for(int i = 0; i < n; i++) {
		L[i].first  = scores[i];
		L[i].second = labels[i];
	}
	std::sort   (L.begin(),L.end());
	std::reverse(L.begin(),L.end());

  	/* Count number of positive and negative examples first */
	int N=0,P=0;
	for(int i = 0; i < n ; i++) {
		if(labels[i] == posclass) P++;
		else N++;
	}
    if( N == 0 || P == 0 )
         mexErrMsgTxt("I only found class 1 in the labels vector ...\n");

    /* Then calculate the actual are under the ROC curve */
	double              A       = 0;
	double              fprev   = INT_MIN; //-infinity
	unsigned long long	FP      = 0, 
                        TP      = 0,
                        FPprev  = 0, 
                        TPprev  = 0;
    
	for(int i = 0 ; i < n; i++) {
		double fi   = L[i].first;
		double label= L[i].second;		
		if(fi != fprev) {
            /* Divide area here already : a bit slower, but gains in precision and avoids overflows */
			A       = A + (trapezoidArea(FP*1.0/N,FPprev*1.0/N,TP*1.0/P,TPprev*1.0/P));
			fprev	= fi;
			FPprev	= FP;
			TPprev	= TP;
		}
		if(label  == posclass)
			TP = TP + 1;
		else
			FP = FP + 1;
	}
    /* Not the same as Fawcett's original pseudocode though I think his contains a typo */
	A = A + trapezoidArea(1.0,FPprev*1.0/N,1.0,TPprev*1.0/P); 
	return A;
}
/* Caculate the trapezoidal area bound by the quad (X1,X2,Y1,Y2)*/
double trapezoidArea(double X1, double X2, double Y1, double Y2) {
	double base   = std::abs(X1-X2);
	double height =     (Y1+Y2)/2.0;
	return (base * height);
}