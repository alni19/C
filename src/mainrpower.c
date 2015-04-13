#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <unistd.h>
#include <string.h>
#include "utilities.h"
#include "power.h"

void *PWR_wrapper(void *pvoidedbag);

int main(int argc, char *argv[])
{
    int code = 0, i, j, initialruns, activeworkers, scheduledjobs;
    int numassets;
    int numdays;
    double lambda = 1;
    double* lb, *ub;
    int r=2;/*factor number*/
    double** matrix;
    double*** matrix_array;
    double* v;
    powerbag **ppbag = NULL, *pbag;
    double scale = 1.0;	
    int quantity = 1, numworkers = 1, theworker;
    char gotone;
    pthread_t *pthethread1;
    pthread_mutex_t output;
    pthread_mutex_t *psynchro_array;
    /*read numassets and numdays*/
    char buffer[100];
    FILE* in = fopen(argv[2],"r");
    fscanf(in,"%s", buffer);
    fscanf(in,"%s", buffer);
    numassets = atoi(buffer);
    fscanf(in,"%s", buffer);
    fscanf(in,"%s", buffer);
    numdays = atoi(buffer);
    fscanf(in,"%s", buffer);
    fscanf(in,"%s", buffer);
    lambda = atof(buffer);
    
    lb = (double*)calloc(numassets,sizeof(double));
    ub = (double*)calloc(numassets,sizeof(double));
    for(i=0;i<numassets;i++){lb[i]=-10000;ub[i]=10000;}

    matrix = (double**) calloc(numassets,sizeof(double*));
    for (i=0; i<numassets; i++) matrix[i] = (double*) calloc(numdays,sizeof(double));/*initialize the time series price matrix*/

    v = (double*)calloc(numassets,sizeof(double));

    /*reads data from CSV file into a matrix*/
    csvread(argv[1],matrix);
    /*
    "/home/alan/gittest/data/dump2.csv"
    */
    calculate_v(matrix,v,numassets,numdays);

  if(argc < 3){
    printf(" usage: rpower data_filename config_filename [-s scale] [-q quantity] [-w workernumber]\n");
    code = 1; goto BACK;
  }

  for(j = 3; j < argc; j++){
    if (0 == strcmp(argv[j],"-s")){
      j += 1;
      scale = atof(argv[j]);
    }
    else if (0 == strcmp(argv[j],"-q")){
      j += 1;
      quantity = atoi(argv[j]);
    }
    else if (0 == strcmp(argv[j],"-w")){
      j += 1;
      numworkers = atoi(argv[j]);
    }
    else{
      printf("bad option %s\n", argv[j]); code = 1; goto BACK;
    }
  }

  printf("will use scale %g and quantity %d: %d workers\n", scale, quantity,
	 numworkers);

  matrix_array = (double***)calloc(quantity,sizeof(double**));

  if( numworkers > quantity ){
    numworkers = quantity; printf(" --> reset workers to %d\n", numworkers);  
  }

  pthread_mutex_init(&output, NULL); /** common to everybody **/

  psynchro_array = (pthread_mutex_t *)calloc(numworkers, sizeof(pthread_mutex_t));
  if(!psynchro_array){
    printf("could not create mutex array\n"); code = NOMEMORY; goto BACK;
  }

  for(j = 0; j < numworkers; j++)
    pthread_mutex_init(&psynchro_array[j], NULL);
  printf("%d",1);
  ppbag = (powerbag **)calloc(numworkers, sizeof(powerbag *));
  if(!ppbag){
    printf("could not create bag array\n"); code = NOMEMORY; goto BACK;
  }

  for(j=0; j<numworkers; j++)
  {
    ppbag[j]=(powerbag*)calloc(1,sizeof(powerbag));
  }
  pthethread1 = (pthread_t *)calloc(numworkers, sizeof(pthread_t));
  if(!pthethread1){
    printf("could not create thread array\n"); code = NOMEMORY; goto BACK;
  }
  /*******till here everything fine*******/
  


  for(j = 0; j < numworkers; j++){
    pbag = ppbag[j];
    pbag->psynchro = &psynchro_array[j];
    pbag->poutputmutex = &output;
    pbag->command = STANDBY;
    pbag->status = PREANYTHING;
    pbag->ID = j;
    pbag->numberAssets = numassets;
    pbag->t = numdays;
    pbag->lambda = lambda;
    pbag->ub = ub;
    pbag->lb = lb;
    pbag->r = r;
    pbag->x=(double*)calloc(numassets,sizeof(double));

    printf("about to launch thread for worker %d\n", j);

    pthread_create(&pthethread1[j], NULL, &PWR_wrapper, (void *) pbag);
  }

  initialruns = numworkers;
  if (initialruns > quantity) initialruns = quantity;

  for(theworker = 0; theworker < initialruns; theworker++){
    pbag = ppbag[theworker];
    matrix_array[theworker] = (double**)calloc(numassets,sizeof(double*));
    for (i=0; i<numassets; i++) matrix_array[theworker][i] = (double*) calloc(numdays,sizeof(double));
    perturb(matrix,matrix_array[theworker],numassets,numdays,v,scale);

    pthread_mutex_lock(&output);
    printf("*****master:  worker %d will run experiment %d\n", theworker, theworker);
    pthread_mutex_unlock(&output);

    /** tell the worker to work **/
    pthread_mutex_lock(&psynchro_array[theworker]);
    pbag->matrix = matrix_array[theworker];
    pbag->command = WORK;
    pbag->status = WORKING;
    pbag->jobnumber = theworker;
    pthread_mutex_unlock(&psynchro_array[theworker]);
  }
  scheduledjobs = activeworkers = initialruns;

  while(activeworkers > 0){
    /** check the workers' status **/
    gotone = 0;
    for(theworker = 0; theworker < numworkers; theworker++){

      pthread_mutex_lock(&psynchro_array[theworker]);
      pbag = ppbag[theworker];
      if(pbag->status == DONEWITHWORK){

	pthread_mutex_lock(&output);
	printf("master:  worker %d is done with job %d\n", pbag->ID, pbag->jobnumber);
	pthread_mutex_unlock(&output);
	if(scheduledjobs >= quantity){
	  /** tell worker to quit **/
	  pthread_mutex_lock(&output);
	  printf("master: telling worker %d to quit\n", pbag->ID);
	  pthread_mutex_unlock(&output);
	  pbag->command = QUIT;
	  pbag->status = QUIT;
	  --activeworkers;
	}
	else {
	  gotone = 1;
	}
      }
      else if(pbag->status == PREANYTHING){
	pthread_mutex_lock(&output);
	printf("master:  worker %d is available\n", theworker);
	pthread_mutex_unlock(&output);
	gotone = 1;
      }

      pthread_mutex_unlock(&psynchro_array[theworker]);
      if(gotone) break;
      usleep(100000);

    }
    /** at this point we have run through all workers **/

    if(gotone){
    /** if we are here, "theworker" can work **/
        pbag = ppbag[theworker];
        matrix_array[scheduledjobs] = (double**)calloc(numassets,sizeof(double*));
        for (i=0; i<numassets; i++) matrix_array[scheduledjobs][i] = (double*) calloc(numdays,sizeof(double));
        perturb(matrix,matrix_array[theworker],numassets,numdays,v,scale);

      pthread_mutex_lock(&output);
      printf("master:  worker %d will run experiment %d\n", theworker, scheduledjobs);
      pthread_mutex_unlock(&output);


      /** tell the worker to work **/
      pthread_mutex_lock(&psynchro_array[theworker]);
      pbag->command = WORK;
      pbag->status = WORKING;
      pbag->jobnumber = scheduledjobs;
      pthread_mutex_unlock(&psynchro_array[theworker]);

      ++scheduledjobs;
    }
  }



  /*  pthread_mutex_lock(&psynchro_array[theworker]);
  pbag->command = QUIT;
  pthread_mutex_unlock(&psynchro_array[theworker]);*/

  pthread_mutex_lock(&output);
  printf("master:  done with loop\n");
  pthread_mutex_unlock(&output);

  /** actually this is bad -- should wait for the threads to be done --
      but how **/
  for(j = 0; j < numworkers; j++){
    pthread_join(pthethread1[j], NULL);
    pthread_mutex_lock(&output);
    printf("master: joined with thread %d\n", j);
    pthread_mutex_unlock(&output);
    pbag = ppbag[j];
    PWRfreespace(&pbag);
  }
  free(ppbag);  


BACK:
  return code;
}


void *PWR_wrapper(void *pvoidedbag)
{
  powerbag *pbag = (powerbag *) pvoidedbag;

  CALLWORKER(pbag);

  return (void *) &pbag->ID;
}





















