#include "glmnet.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include <vector>
using namespace std;

float RandDbl(float d){
	//random float btw 0 and d
	return (rand()*d)/(RAND_MAX+1.f);
}

float NormRand(float &next_gaussian){
	float fac, rsq, v1, v2;

	do {
		v1 = 2.0*RandDbl(1)-1.0;
		v2 = 2.0*RandDbl(1)-1.0;
		rsq = v1*v1+v2*v2;
	} while (rsq >= 1.0 || rsq == 0.0);
	fac = sqrt(-2.0*log(rsq)/rsq);
	next_gaussian=v1*fac;
	return v2*fac;
}

void GetObsLogRegN(float *x, float &y, int M, int kstar, float alpha){
	float a2=alpha/2;
	if (RandDbl(1)<0.5f)
		y=0;
	else
		y=1;
	for (int m=0;m<kstar;m+=2){
		if (y==1){
			x[m]=(NormRand(x[m+1])+a2)/(5+a2);
			x[m+1]=(x[m+1]+a2)/(5+a2);
		}
		else{
			x[m]=(NormRand(x[m+1])-a2)/(5+a2);
			x[m+1]=(x[m+1]-a2)/(5+a2);
		}
		if(x[m]>1) x[m]=1; 
		if (x[m]<-1) x[m]=-1;
		if(x[m+1]>1) x[m+1]=1; 
		if (x[m+1]<-1) x[m+1]=-1;
	}

	for (int m=kstar;m<M;m+=2){
		x[m]=NormRand(x[m+1])*0.2f;
		x[m+1]*=0.2f;
		if(x[m]>1) x[m]=1; 
		if (x[m]<-1) x[m]=-1;
		if(x[m+1]>1) x[m+1]=1; 
		if (x[m+1]<-1) x[m+1]=-1;
	}
}

int TestGlmnet(int N, int M,int kstar, float alpha, int seed){
	vector<float> x,xall,y,beta;
	vector<vector<float> > allbeta;
	float sumb;
	int i,j,m,nnz,n1;
	GlmnetResult rez;
	sumb=2;
	x.resize(M);
	xall.resize(M*N);
	y.resize(N);
	srand(seed);
	for (i=0;i<N;i++){
		GetObsLogRegN(&x[0],y[i],M,kstar,alpha);
		for (j=0;j<M;j++)
			xall[j*N+i]=x[j];
	}
	DoGlmnet(rez,&xall[0],&y[0],N,M,2);
	rez.GetBetas(allbeta,0);
	for (i=0;i<rez.nLambdas;i++)
		if (rez.nNonZeroCoeff[i]==kstar){
			beta=allbeta[i];
			break;
		}
	if (beta.empty())
		return 0;
	n1=0;
	for (m=0;m<kstar;m++)
		if (beta[m]!=0)
			n1++;
	nnz=0;
	for (m=kstar;m<M;m++)
		if (beta[m]!=0)
			nnz++;
	if (n1==kstar&&nnz==0)
		return 1;
	return 0;
}

int main(int argc, char* argv[]){
	if (argc<=3){
		printf("Usage: testglmnet nObs nVar Alpha [kstar] \n");
		return 1;
	}

	int N=atoi(argv[1]);
	int M=atoi(argv[2]);
	float alpha=atof(argv[3]);
	vector<int> result;
	char erase[7];
	int r,nRuns=200,sum,kstar=10;
    memset(erase,'\b',sizeof(erase));

	if (argc>4)
		kstar=atoi(argv[4]);

	printf("N=%d,M=%d,kstar=%d,alpha=%1.2f,nRuns=%d\n",N,M,kstar,alpha,nRuns);
	result.resize(nRuns);
	int nthreads=omp_get_max_threads();
//#pragma omp parallel for //num_threads(4)
	for (r=0;r<nRuns;r++){
        int thread_id = omp_get_thread_num();
		time_t start;
		time(&start);
		result[r]=TestGlmnet(N,M,kstar,alpha,r);
		//printf("%d",result[r]);
		if (r > 0)
            fwrite(erase,1,sizeof(erase),stdout);
		printf("%d  (% 3ds)",result[r],time(NULL)-start);
	}
	sum=0;
	for (r=0;r<nRuns;r++) 
		sum+=result[r];
	char st[255];
	sprintf(st,"glmnet_%dx%d_%1.2f.txt",N,M,alpha);
	FILE *f=fopen(st,"a+");
	if (f!=0){
		for (r=0;r<nRuns;r++){
			if (r<nRuns-1)
				fprintf(f,"%d ",result[r]);
			else
				fprintf(f,"%d\n",result[r]);
		}
		fclose(f);
	}
	printf("\ndet=%d,perc=%1.2f\n",sum,sum/(float)nRuns*100);
	return 0;
}
