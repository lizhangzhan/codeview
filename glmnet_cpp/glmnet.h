#ifndef GLMNET_H
#define GLMNET_H

#include<vector>

#ifdef __cplusplus
extern "C" {
#endif
//c   x, ix, jx = predictor data matrix in compressed sparse column format
//c   ka = algorithm flag
//c      ka=1 => covariance updating algorithm
//c      ka=2 => naive algorithm
//c   parm = family member index (0 <= parm <= 1)
//c        = 0.0 => ridge
//c        = 1.0 => lasso
//c   no = number of observations
//c   ni = number of predictor variables
//c   y(no) = response vector
//c   w(no)= observation weights
//c   jd(jd(1)+1) = predictor variable deletion flag
//c      jd(1) = 0  => use all variables
//c      jd(1) != 0 => do not use variables jd(2)...jd(jd(1)+1)
//c   vp(ni) = relative penalties for each predictor variable
//c      vp(j) = 0 => jth variable unpenalized
//c   ne = maximum number of variables allowed to enter largest model
//c        (stopping criterion)
//c   nx = maximum number of variables allowed to enter all models
//c        along path (memory allocation, nx > ne).
//c   nlam = (maximum) number of lamda values
//c   flmin = user control of lamda values (>=0)
//c      flmin < 1.0 => minimum lamda = flmin*(largest lamda value)
//c      flmin >= 1.0 => use supplied lamda values (see below)
//c   ulam(nlam) = user supplied lamda values (ignored if flmin < 1.0)
//c   thr = convergence threshold for each lamda solution.
//c      iterations stop when the maximum standardized coefficient
//c      change from the previous iteration is less than thr
//c      (suggested value, thr=1.0e-4)
//c   isd = standarization flag:
//c      isd = 0 => regression on original predictor variables
//c      isd = 1 => regression on standardized predictor variables
//c      Note: output solutions always reference original
//c            variables locations and scales.
//c
//c output:
//c
//c   lmu = actual number of lamda values (solutions)
//c   a0(lmu) = intercept values for each solution
//c   ca(nx,lmu) = compressed coefficient values for each solution
//c   ia(nx) = pointers to compressed coefficients
//c   nin(lmu) = number of compressed coefficients for each solution
//c   rsq(lmu) = R**2 values for each solution
//c   alm(lmu) = lamda values corresponding to each solution
//c   nlp = total passes over the data summed over all lamda values
//c   jerr = error flag:
//c      jerr  = 0 => no error
//c      jerr != 0 => fatal error - no output returned
//c         jerr < 7777 => memory allocation error
//c         jerr = 7777 => all used predictors have zero variance
//c         jerr = 10000 => maxval(vp) <= 0.0
void ELNET(int *ka, float *parm, int *no, int *ni, float *x, float *y, float *w, int *jd, 
		   float *vp, int *ne, int *nx, int *nlam, float *flmin, float *ulam, float *thr, int *isd, int *lmu,
		   float *a0, float *ca, int *ia, int *nin, float *rsq, float *alm, int *nlp, int *jerr);

//c   x, ix, jx = predictor data matrix in compressed sparse column format
//c   parm, no, ni, jd, vp, ne, nx, nlam, flmin, ulam, thr, isd, same as above.
//c
//c   nc = number of classes (distinct outcome values)
//c        nc=1 => binomial two-class logistic regression
//c            (all output references class 1)
//c   y(no,max(2,nc)) = number of each class at each design point(overwritten)
//c   maxit = maximum number of iterations allowed for any lamda value
//c           (suggested value, maxit = 100)
//c   kopt = optimization flag
//c      kopt = 0 => Newton-Raphson
//c      kpot = 1 => modified Newton-Raphson (recommended)
//c
//c
//c output:
//c
//c   lmu, ia, nin, alm, nlp, same as above
//c
//c   a0(nc,lmu) = intercept values for each class at each solution
//c   ca(nx,nc,lmu) = compressed coefficient values for each class at
//c                each solution
//c   dev(lmu) = fraction of explained devience for each solution
//c   jerr = error flag
//c      jerr = 0  => no error
//c      jerr > 0 => fatal error - no output returned
//c         jerr < 7777 => memory allocation error
//c         jerr = 7777 => all used predictors have zero variance
//c         jerr = 8000 + k => null probability < 1.0e-5 for class k
//c         jerr = 9000 + k => null probability for class k
//c                            > 1.0 - 1.0e-5
//c         jerr = 10000 => maxval(vp) <= 0.0
//C      jerr < 0 => non fatal error - partial output returned
//c         jerr = -k => convergence for kth lamda value not reached
//c            after maxit (see above) iterations. Solutions for
//c            larger lamdas returned
//c         jerr = -10000-k => number of non zero coefficients along path
//c            exceeds nx (see above) at kth lamda value. Solutions for
//c            larger lamdas returned
void LOGNET(float *parm, int *no, int *ni, int *nc, float *x, float *y, int *jd, float *vp, int *ne,
			int *nx, int *nlam, float *flmin, float *ulam, float *thr, int *isd, int *maxit, int *kopt, int *lmu,
			float *a0, float *ca, int *ia, int *nin, float *dev, float *alm, int *nlp, int *jerr);

#ifdef __cplusplus
}
#endif
struct GlmnetResult{
	int nLambdas,nlp,jerr,nVar,nClasses,nx;
	std::vector<float> intercepts,coeffs,alm,dev;
	std::vector<int> ia,nNonZeroCoeff;
	void GetBeta(std::vector<float> &beta, int idxLambda, int cls);
	void GetBetas(std::vector<std::vector<float> > &beta, int cls);
};

int DoGlmnet(GlmnetResult &out, float *x, float *y, int nObs, int nVar, int nClasses);

#endif
