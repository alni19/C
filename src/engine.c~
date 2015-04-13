
#include <stdio.h>
#include <stdlib.h>
#include <string.h> 
#include <gurobi_c.h>
#include <pthread.h>
#include <unistd.h>
#include "utilities.h"
#include "power.h"

int engine(int numassets, int numfactors, 
	     double *ub, double *lb, double *mu, double *sigma2, 
	     double *V, double *F,double lambda,pthread_mutex_t *poutputmutex,int ID, powerbag* pbag)
{
  int retcode = 0;
  GRBenv   *env = NULL;
  GRBmodel *model = NULL;
  int n, i, j, k;
  double *x;
  int *qrow, *qcol, Nq;
  double *qval;
  double *obj;
  int *cind;
  double rhs;
  char sense;
  double *cval;
  int numnonz;
  double sharperatio;

  char **names, bigname[100];

  pthread_mutex_lock(poutputmutex);
  printf("ID %d running solver engine\n",ID);
  pthread_mutex_unlock(poutputmutex);

  n = numassets + numfactors;

  retcode = GRBloadenv(&env, "factors.log");
  if (retcode) goto BACK;

 /* Create initial model */
  retcode = GRBnewmodel(env, &model, "factor", n, 
                      NULL, NULL, NULL, NULL, NULL);
  if (retcode) goto BACK;

  names = (char **) calloc(n, sizeof(char *));

  /** next we create the remaining attributes for the n columns **/
  x     = (double *) calloc (n, sizeof(double));
  obj = (double* ) calloc (1,sizeof(double));

  for(j = 0; j < numassets; j++){
    names[j] = (char *)calloc(3, sizeof(char));
    if(names[j] == NULL){
      retcode = 1; goto BACK;
    }
    sprintf(names[j],"x%d", j);
  }
  for(j = numassets; j < n; j++){
    names[j] = (char *)calloc(3, sizeof(char));
    if(names[j] == NULL){
		  retcode = 1; goto BACK;
    }
    sprintf(names[j],"y%d", j - numassets);
  }
  /* initialize variables */
  for(j = 0; j < n; j++){
    retcode = GRBsetstrattrelement(model, "VarName", j, names[j]);
    if (retcode) goto BACK;

    retcode = GRBsetdblattrelement(model, "Obj", j, -mu[j]);
    if (retcode) goto BACK;

    retcode = GRBsetdblattrelement(model, "LB", j, lb[j]);
    if (retcode) goto BACK;

    retcode = GRBsetdblattrelement(model, "UB", j, ub[j]);
    if (retcode) goto BACK;
  }

  /** next, the quadratic -- there are numassets + numfactors*numfactors nonzeroes: 
      numassets residual variances plus the numfactors x numfactors
      factor covariance matrix**/

  Nq = numassets + numfactors*numfactors;
  qrow = (int *) calloc(Nq, sizeof(int));  /** row indices **/
  qcol = (int *) calloc(Nq, sizeof(int));  /** column indices **/
  qval = (double *) calloc(Nq, sizeof(double));  /** values **/

  if( ( qrow == NULL) || ( qcol == NULL) || (qval == NULL) ){
    pthread_mutex_lock(poutputmutex);
    printf("could not create quadratic\n");
    pthread_mutex_unlock(poutputmutex);
    retcode = 1; goto BACK;

  }

  for (j = 0; j < numassets; j++){
    qval[j] = lambda*sigma2[j];
    qrow[j] = qcol[j] = j;
  }
  for (i = 0; i < numfactors; i++){
    for (j = 0; j < numfactors; j++){
      k = i*numfactors + j;
      qval[k + numassets] = lambda*F[k];
      qrow[k + numassets] = numassets + i;
      qcol[k + numassets] = numassets + j;
    }
  }
  retcode = GRBaddqpterms(model, Nq, qrow, qcol, qval);
  if (retcode) goto BACK;

  /** now we will add one constraint at a time **/
  /** we need to have a couple of auxiliary arrays **/

  cind = (int *)calloc(n, sizeof(int));  /** n is over the top since no constraint is totally dense;		     but it's not too bad here **/
  cval= (double *)calloc(n, sizeof(double));
  if(!cval){
    pthread_mutex_lock(poutputmutex);
    printf("cannot allocate cval\n"); retcode = 2; goto BACK;
    pthread_mutex_unlock(poutputmutex);
  }
  for(i = 0; i < numfactors; i++){
    for(j = 0; j < numassets; j++){
      cval[j] = V[i*numassets + j];
      cind[j] = j;
    }
    cind[numassets] = /* j */ numassets + i;
    cval[numassets] = -1;
    numnonz = numassets + 1;
    rhs = 0;
    sense = GRB_EQUAL;

    sprintf(bigname,"factor%d",i);
    retcode = GRBaddconstr(model, numnonz, cind, cval, sense, rhs, bigname);
    if (retcode) goto BACK;

  }

  retcode = GRBupdatemodel(model);
  if (retcode) goto BACK;

  /** optional: write the problem **/

  retcode = GRBwrite(model, "factorqp.lp");
  if (retcode) goto BACK;

  retcode = GRBoptimize(model);
  if (retcode) goto BACK;

  /** get solution **/

  retcode = GRBgetdblattrarray(model, GRB_DBL_ATTR_X, 0, n, x);

  for(j = 0; j < numassets; j++){
    (pbag->x)[j]=x[j];
  }
  if(retcode) goto BACK;
  retcode = GRBgetdblattr(model, GRB_DBL_ATTR_OBJVAL, obj);
  pbag->obj=*obj;
  if(retcode) goto BACK;  
  sharperatio = sharpe_ratio(pbag->matrix, pbag->x, numassets, pbag->t);

  /** now let's see the values and analysis**/
  pthread_mutex_lock(poutputmutex);
  for(j = 0; j < n; j++){
    printf("%s = %g\n", names[j], x[j]);
  }
  printf("\n*****************Analysis****************\n");
  printf("Objective = %g\n",*obj);
  printf("Portfolio's sharpe ratio = %g\n", sharperatio);
  printf("*****************************************\n\n");
  pthread_mutex_unlock(poutputmutex);
  
  GRBfreemodel(model);
  GRBfreeenv(env);
  
  



 BACK:
  pthread_mutex_lock(poutputmutex);
  printf("ID %d engine exits with code %d\n", ID, retcode);
  pthread_mutex_unlock(poutputmutex);
  return retcode;
}
