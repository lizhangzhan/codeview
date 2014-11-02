#include "glmnet.h"
#include<vector>

using namespace std;

int DoGlmnet(GlmnetResult &out, float *x, float *y, int nObs, int nVar, int nClasses){
	//  each obs in x is on a column
	float parm,flmin,thr;
	int nx,nlam,isd,maxit,kopt,yi,nMaxVar;
	int jd[2];
	vector<float> vp,y1,ulam;

	//    Output
	//int jerr;
	//int lmu,nlp,jerr;
	//vector<float> a0,ca,alm,dev;
	//vector<int> ia,nin;

	y1.assign(nObs*nClasses,0);
	for (int i=0;i<nObs;i++){
		yi=(int)y[i];
		y1[yi*nObs+i]=1;
	}
	parm=1.f;// l1 penalty
    jd[0] = jd[1]= 0;
	nMaxVar=nVar;
	if (nMaxVar==0)
		nMaxVar=nVar+1;
	vp.assign(nVar,1); // individual lambdas
	nx=min((int)(nMaxVar*1.2),nVar); // max number of nonzero betas overall
    nlam = 100; // number of lambdas
	if (nObs<nVar)
		flmin = 0.05f;
	else
		flmin = 0.0001f; // smallest lambda

    ulam.assign(2,0); 
    thr = 0.0001f;
	isd = 1;
	kopt=1;

    maxit = 100;

	out.nVar=nVar;
	out.nx=nx;
	out.nClasses=nClasses;
    out.nLambdas = 0;
    out.nlp = 0;
    out.jerr = 0;
	out.ia.assign(nx,0);
	out.nNonZeroCoeff.assign(nlam,0);
	out.alm.assign(nlam,0);
	out.intercepts.assign(nClasses*nlam,0);
	out.coeffs.assign(nx*nClasses*nlam,0);
	out.dev.assign(nlam,0);
	
	LOGNET(&parm,&nObs,&nVar,&nClasses,x,&y1[0],jd,&vp[0],&nMaxVar,&nx,&nlam,&flmin,&ulam[0],&thr,&isd,&maxit,&kopt,
		&out.nLambdas,&out.intercepts[0],&out.coeffs[0],&out.ia[0],&out.nNonZeroCoeff[0],&out.dev[0],&out.alm[0],&out.nlp,&out.jerr);
	return out.jerr;
}

void GlmnetResult::GetBeta(std::vector<float> &beta, int idxLambda, int cls){
	beta.resize(nVar);
	for (int j=0;j<nVar;j++)
		beta[j]=coeffs[idxLambda*nClasses*nx+j+cls];
}

void GlmnetResult::GetBetas(std::vector<std::vector<float> > &beta, int c){
	beta.resize(nLambdas);
	for (int i=0;i<nLambdas;i++)
		GetBeta(beta[i],i,c);
}
