#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>
#include "utilities.h"



char does_it_exist(char *filename)
{
   struct stat buf;

/* the function stat returns 0 if the file exists*/

   if (0 == stat(filename, &buf)){
      return 1;
   }
   else return 0;
}

void gotosleep(int numseconds)
{
   sleep(numseconds);
}

void erasefile(char *filename)
{
  remove(filename);
}


double drawnormal(void)
{
  double U1, U2, drawn, pi;

  pi = 3.141592653589793;

  U1 = rand()/((double) RAND_MAX);
  U2 = rand()/((double) RAND_MAX);

  drawn = sqrt(-2*log(U1))*cos(2*pi*U2);

  return drawn;
}

void normal_sum0_init(double* arr, int N)
{
    srand(time(NULL));
    int i;
    double sum=0,average;
    for (i=0;i<N;i++)
    {
        arr[i]=drawnormal();
        sum=sum+arr[i];
    }
    average=sum/N;
    for (i=0;i<N;i++)
    {
        arr[i]=arr[i]-average;
    }
}

void split(double* price,char* str){
    double num;
    int i=0;
    char * pch;
    pch = strtok (str," ,-");
    while (pch != NULL)
    {
        price[i]=atof(pch);
        //printf("%f\n",price[i]);
        i++;
        pch = strtok (NULL, " ,-");
    }
}
void calculate_v(const double** matrix, double* v,int numassets, int numdays)
{
	int i,j;
	double mu;
	double* r;
	r=calloc(numdays-1,sizeof(double));
	for (i=0;i<numassets;i++)
	{
        mu=0;
		for (j=0;j<numdays-1;j++)
		{
            r[j]=(matrix[i][j+1]-matrix[i][j])/matrix[i][j];
            mu+=r[j];
        }
        mu=mu/(numdays-1);
        v[i]=0;
        for (j=0;j<numdays-1;j++)
            v[i]+=(r[j]-mu)*(r[j]-mu);
        v[i]=sqrt(v[i]/(numdays-1));
    }
}

void perturb(const double** matrix, double** target, int numassets, int numdays,const double* v)
{
    int i,j;
    double* eps;
    eps = (double*)calloc(numdays-1,sizeof(double));
    normal_sum0_init(eps,numdays-1);
    for(i=0;i<numassets;i++)
    {
        target[i][0]=matrix[i][0];
        for (j=0;j<numdays-1;j++)
        {
            target[i][j+1]=target[i][j]*(1+(matrix[i][j+1]-matrix[i][j])/matrix[i][j]+v[i]*eps[j]);
        }
    }
}
