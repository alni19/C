#include <pthread.h>
#include <unistd.h>
#include "utilities.h"
#include "power.h"
int engine(int numassets, int numfactors, 
	     double *ub, double *lb, double *mu, double *sigma2, 
	     double *V, double *F,double lambda,pthread_mutex_t *poutputmutex );

void poweralg(int ID, int n, double *vector, double *newvector, double *matrix, powerbag *pbag, double *eigenVal)
{
  int k;
  double error, tolerance;
  tolerance = .000001/n;

  for(k = 0; ; k++)
  {
    poweriteration(ID,k, n, vector, newvector, matrix ,&error, pbag->poutputmutex, eigenVal);
    if(error < tolerance)
    {
	pthread_mutex_lock(pbag->poutputmutex);
	printf(" ID %d converged to tolerance %g! on job %d\n", ID, tolerance,
	       pbag->jobnumber); 
	printf(" ID %d top eigenvalue  %g!\n", ID, *eigenVal);
	pthread_mutex_unlock(pbag->poutputmutex);
        break;
    }
  }
}


void poweriteration(int ID, int k, int n, double *vector, double *newvector, double *matrix, double *perror,pthread_mutex_t *poutputmutex, double *eigVal)
{
  double norm2 = 0, mult, error;
  int i, j;
  for(i = 0; i <n; i++){
    newvector[i] = 0;
    for (j = 0; j < n; j++){
      newvector[i] += vector[j]*matrix[i*n + j];
    }
  }

  norm2 = 0;
  for(j = 0; j < n; j++) norm2 += newvector[j]*newvector[j];

  mult = 1/sqrt(norm2);

  for(j = 0; j < n; j++) newvector[j] = newvector[j]*mult;

  PWRcompute_error(n, &error, newvector, vector);

  if(0 == k%100)
  {
    pthread_mutex_lock(poutputmutex);
    printf("ID %d at iteration %d, norm is %g, ", ID, k, 1/mult);
    printf("  L1(error) = %.9e\n", error);
    pthread_mutex_unlock(poutputmutex);
  }
  
  *eigVal = 1/mult;
  *perror = error;

  /** will need to map newvector into vector if not terminated **/

  for(j = 0; j < n; j++) vector[j] = newvector[j];

}

void PWRcompute_error(int n, double *perror, double *newvector, double *vector)
{
  int j;
  double error;

  error = 0;

  for(j = 0; j < n; j++){
    error += fabs(newvector[j] - vector[j]);
  }

  *perror = error;

}

void newOmegaVector(int n, double *omegaVector, double *vector)
{
  double sum =0;
  int i;
  for(i=0;i<n;i++)
  {
    sum += omegaVector[i]*vector[i];  
  }
  for(i=0;i<n;i++)
  {
    omegaVector[i] = omegaVector[i]-sum*vector[i];
  } 
}


void newMatrix(int n,double *matrix, double *vector,double eigenVal)
{ int i,j;
  for(i=0;i<n;i++)
  {
	 for (j=0;j<n;j++)
     {
		 matrix[i*n+j]=matrix[i*n+j]-eigenVal*vector[i]*vector[j];
  	 }
  }
}

void PWRfreespace(powerbag **ppbag)
{
  powerbag *pbag = *ppbag;
  if(pbag == NULL) return;

  PWRfree3(&pbag->matrix, pbag->numberAssets);
  PWRfree(&pbag->ub);
  PWRfree(&pbag->lb);

  free(pbag);
  *ppbag = NULL;
}

void PWRfree(double **ppaddress)
{
  double *paddress = *ppaddress;
  if(paddress == NULL) return;
  printf("freeing double array at %p\n", (void *) paddress);
  free(paddress);
  *ppaddress = NULL; /** prevents double freeing **/
}
void PWRfree3(double ***pmatrix, int n)
{
   int i;
   for(i=0;i<n;i++)
   free((*pmatrix)[i]);
   free(*pmatrix);
   *pmatrix = NULL;
}


/********************************core function**********************************/

void CALLWORKER(powerbag *pbag)
{
  /***********Part Zero: Initialization***********/
  int ID;
  int retcode = 0;
  int i,t,j,n,r;
  int waitcount;
  char letsgo = 0;
  double *matrix = NULL, *vector = NULL, *newvector = NULL;
  double *eigenVal=NULL;
  double **eigenVector=NULL;
  double *omegaVector=NULL;
  double *Qmatrix=NULL;
  double *Vtranspose = NULL;
  double *F =NULL;
  double *residual= NULL;
  double *mu = NULL;
  double **timeSeries = NULL;
  double *ub = NULL;
  double *lb = NULL;
  double lambda;  
ub = pbag->ub;
lb = pbag->lb;  
lambda = pbag->lambda;
  
  timeSeries=pbag->matrix;
  r=pbag->r;
  n=pbag->numberAssets;
  t = pbag->t;
  ID = pbag->ID;
  pbag->status = WAITING;
  
  pthread_mutex_lock(pbag->poutputmutex);
  printf("ID %d starts\n", pbag->ID);
  pthread_mutex_unlock(pbag->poutputmutex);

  for(;;)
  {
    pthread_mutex_lock(pbag->poutputmutex);
    printf(" ID %d in big loop\n", pbag->ID);
    pthread_mutex_unlock(pbag->poutputmutex);

    letsgo = 0;
    waitcount = 0;
    while(letsgo == 0)
    {
      /** wait until WORK signal **/
      sleep(1);

      pthread_mutex_lock(pbag->psynchro);
      if(pbag->command == WORK)
      {
	letsgo = 1;
      }
      else if(pbag->command == QUIT)
	letsgo = 2;
      pthread_mutex_unlock(pbag->psynchro);

      if (letsgo == 2) 
	goto DONE;

      if(0 == waitcount%2)
      {
	pthread_mutex_lock(pbag->poutputmutex);
	printf("ID %d: wait %d for signal to start working\n", pbag->ID, waitcount);
	pthread_mutex_unlock(pbag->poutputmutex);

      }
      ++waitcount;

    }
}
    pthread_mutex_lock(pbag->poutputmutex);
    printf("ID %d: got signal to start working\n", pbag->ID);
    pthread_mutex_unlock(pbag->poutputmutex);

  /*allocate eigen value array, omegavector array, eigenvector matrix*/
  matrix = (double*) calloc(n*n,sizeof(double));
  vector = (double*) calloc(n,sizeof(double));
  newvector = (double*) calloc(n,sizeof(double));
  eigenVal = (double*) calloc(r,sizeof(double));
  omegaVector = (double*) calloc(n,sizeof(double));
  eigenVector = (double **) calloc(r,sizeof(double *));
  mu = (double *) calloc(n,sizeof(double));

  
  /********calculate mu and cov as matrix***********/
  calcMuCov(timeSeries,n,t,mu,matrix);  
  
  for(i=0;i<r;i++)
  {
   eigenVector[i] = (double*)calloc(n,sizeof(double));
  }

  Qmatrix=(double*) calloc(n*n,sizeof(double));

  for(i=0;i<n;i++)
  {
	  for(j=0;j<n;j++)
	  {
          Qmatrix[i*n+j]=matrix[i*n+j];
	  }
  }

  Vtranspose = (double*) calloc(n*r,sizeof(double));
  F = (double *) calloc(r*r,sizeof(double));
  residual = (double*)calloc(n,sizeof(double));    

  /*********first part of work:PCA********/  

  /** generate a random vector**/  
  for(j = 0; j < n; j++){ 
    vector[j] = rand()/((double) RAND_MAX);
  }
  
  /*record original omega vector*/
  for (i=0;i<n;i++)
  {
	  omegaVector[i] = vector[i]; 
  }

  /*power method recursion to find eigen decomposition with factor r*/
  for(i=0;i<r;i++)
  {
	  poweralg(ID,n, vector, newvector, matrix,pbag, &(eigenVal[i]));
	  for(j=0;j<n;j++)
	  {
		  eigenVector[i][j] = vector[j];
	  }	  
	  newOmegaVector(n,omegaVector,vector);
          /*let vector=new omega vector*/
	  newMatrix(n,matrix,vector,eigenVal[i]);
          for(j=0;j<n;j++)
          {
              vector[j] = omegaVector[j];
          }
  }
  
  for(i=0;i<n;i++)
  {
	  for(j=0;j<r;j++)
	  {
           Vtranspose[r*i+j]=eigenVector[j][i];
	  }
  }  
  
  for(i=0;i<r;i++)
  {
	  for(j=0;j<r;j++)
	  {
		  if(i==j) F[i*r+j]=eigenVal[i];
                  else F[i*r+j]=0;
		  
	  }
  }

  for(i=0;i<n;i++)
  {
     residual[i] = fabs(Qmatrix[i*n+i]-matrix[i*n+i]);
  }  

  


  /*********Part two: Grubi Optimization*********/

  pthread_mutex_lock(pbag->poutputmutex);
  printf("numassets: %d; numfactors: %d\n", n, r);
  pthread_mutex_unlock(pbag->poutputmutex);
  
  retcode = engine(n,r,ub,lb,mu,residual,Vtranspose,F,lambda,pbag->poutputmutex);
  if(retcode) goto BACK;  
  /**********Done with work**********************/
   pthread_mutex_lock(pbag->psynchro);
   pbag->status = DONEWITHWORK;
   pbag->command = STANDBY;
   pthread_mutex_unlock(pbag->psynchro);
   



  DONE:
  pthread_mutex_lock(pbag->poutputmutex);
  printf(" ID %d quitting\n", pbag->ID);
  pthread_mutex_unlock(pbag->poutputmutex);  

  BACK:
  pthread_mutex_lock(pbag->poutputmutex);
  printf(" retcode= %d\n", retcode);
  pthread_mutex_unlock(pbag->poutputmutex);
}













