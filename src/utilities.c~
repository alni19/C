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
    int i;    
    double sum=0,average;
    srand(time(NULL));
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
    int i=0;
    char * pch;
    pch = strtok (str," ,-");
    while (pch != NULL)
    {
        price[i]=atof(pch);
        i++;
        pch = strtok (NULL, " ,-");
    }
}



void calculate_v(double** matrix, double* v,int numassets, int numdays)
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

void perturb(double** matrix, double** target, int numassets, int numdays,double* v,double scale)
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
    for(i=0;i<numassets;i++)
        for (j=0;j<numdays;j++)
            target[i][j]=(1-scale)*matrix[i][j]+scale*target[i][j];
}

void calcMuCov(double **matrix, int n, int t, double *mu, double *cov)
{
  int i,j,k;
  double sum;
  double *returnMatrix;
  returnMatrix = (double*) calloc(n*(t-1),sizeof(double));



  for(i=0;i<n;i++)
  {
     sum=0;
     for(j=0;j<(t-1);j++)
     {
        returnMatrix[i*(t-1)+j] = (matrix[i][j+1]-matrix[i][j])/(matrix[i][j]); 
        sum = sum + returnMatrix[i*(t-1)+j];
     }
     mu[i] = sum/(t-1);
  }

  for(i=0;i<n;i++)
  {
     for(j=0;j<n;j++)
     { 
        sum=0;
        for(k=0;k<(t-1);k++)
        {
           sum = sum + (returnMatrix[i*(t-1)+k]-mu[i])*(returnMatrix[j*(t-1)+k]-mu[j]);
        } 
        cov[i*n+j]=sum/(t-1);
     }
  }
}


void csvread(const char *filename,double** matrix)
{
    char line[6000];
    int i=0;
    FILE* stream = fopen(filename, "r");
    printf("%s","Loading CSV...");
    while (fgets(line, 6000, stream))
    {
        split(matrix[i],line);
        i++;
    }
    printf("%s","Loading ends.");
}

double sharpe_ratio(double** matrix, double *portfolio, int numassets, int numdays)
{
    int i,j;
    double* newprice = (double* ) calloc(numdays, sizeof(double));
    double sum=0,sum_square=0;
    double sigma, mu;
    for (i=0; i<numdays; i++)
    {
        for (j=0; j<numassets;j++)
            {
                newprice[i]+=matrix[j][i]*portfolio[j];
            }    
    }
    for (i=0; i<numdays-1; i++)
    {
        sum += (newprice[i+1]-newprice[i])/newprice[i];
        sum_square += (newprice[i+1]-newprice[i])*(newprice[i+1]-newprice[i])/newprice[i]/newprice[i];
    }
    mu = sum/(numdays-1);
    sigma = sqrt(sum_square/(numdays-1)-mu*mu);
    return mu/sigma*sqrt(252);
    free(newprice);     
}
