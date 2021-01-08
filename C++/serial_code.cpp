#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
#include <cstdlib>
#include <stdlib.h>
#include <iomanip> 
#include <cmath>
#include <time.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#define NINT(a) ((a) >= 0.0 ? (int)((a)+0.5) : (int)((a)-0.5))
const double PI  =3.141592653589793238463;

using std::cin;
using std::cout;
using std::endl;
using std::vector;

double rand_normal(double mean, double stddev)	//Box muller method
{
    static double n2 = 0.0;
    static int n2_cached = 0;
    if (!n2_cached)
    {
        double x, y, r;
        do
        {
            x = 2.0*rand()/RAND_MAX - 1;
            y = 2.0*rand()/RAND_MAX - 1;

            r = x*x + y*y;
        }
        while (r == 0.0 || r > 1.0);
        {
            double d = sqrt(-2.0*log(r)/r);
            double n1 = x*d;
            n2 = y*d;
            double result = n1*stddev + mean;
            n2_cached = 1;
            return result;
        }
    }
    else
    {
        n2_cached = 0;
        return n2*stddev + mean;
    }
}

struct material
{
	int num = 10;  //number of material 
	float pos[2][100]; //the centre of under 100 inclusions stored
	float diff = 1.0;
};

struct syst
{
	long int nsteps = 200; // lets assume 10 random walk steps
	int n_part = 100;  // here we start with number of particles as 2
	float time = 0.1; // time increment is 3 pico sec
	float dx = 0.1; // scaling factor
	int iter = 3; // number of iteration (parallelised)
	float box[2] = {10.0,10.0}; // lengths of x and y
};

void BC(float& x, float& y, syst param)
{
    //reflective boundary conditions
    if( x < 0.0 )
        x = -x;

    else if( x > param.box[0] )
        x = 2*param.box[0] - x;
  
    //perodic boundary conditions
    if( y  > param.box[1] )
        y = y - param.box[1];

    else if( y < 0.0 )
        y = y + param.box[1];

}

int main(int argc, char** argv) 
{
	clock_t t1,t2;
    t1=clock();
    srand(time(NULL));

    syst param;
    int i,k,ix,iy,it,n_part,i_iter;
    const int rows = param.nsteps;
    const int cols = NINT(param.box[0]/param.dx)+1;
    float hist[rows][cols];
    float hist_avg[rows][cols];
    material matrix; 
    matrix.diff = 1; 
    float sigma = sqrt(matrix.diff*param.time); 

    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
        {
            hist[i][j] = 0.0;
            hist_avg[i][j] = 0.0;
        }

    float *u = (float*) malloc(sizeof(float)*(param.n_part*2*param.nsteps+1)); 
    
    for(int i=0;i<(param.nsteps*param.n_part);i++)
        u[i]=1.0;

    for(int i=(param.nsteps*param.n_part);i<param.nsteps*2*param.n_part;i++)
        u[i]=-1.0;
        
    i_iter = 0;
    while( i_iter < param.iter)
    {
        gsl_rng * r;
        const gsl_rng_type * T;
        T = gsl_rng_gfsr4;
        r = gsl_rng_alloc (T);
        gsl_rng_set (r, (i_iter+1)*rand());

	    float **x = (float**)malloc(sizeof(float*) * 2);
	    for (int i = 0; i < 2; i++) 
	        x[i] = (float *)malloc(sizeof(float) * (param.n_part*2*param.nsteps));

	    for(int i=0;i<(param.nsteps*param.n_part);i++) 
	        x[0][i] = 0.0;

	    for(int i=(param.nsteps*param.n_part);i<2*param.nsteps*param.n_part;i++) 
	        x[0][i] = param.box[0];

	    for(int i=0;i<2*param.nsteps*param.n_part;i++) 
	        x[1][i] = gsl_rng_uniform(r)*param.box[1];

        n_part = 2*param.n_part;

        for(int i=0;i<(param.nsteps*param.n_part);i++) 
            x[0][i] = 0;

        for(int i=(param.nsteps*param.n_part);i<2*param.nsteps*param.n_part;i++) 
            x[0][i] = param.box[0];

        for(int i=0;i<2*param.nsteps*param.n_part;i++) 
            x[1][i] = gsl_rng_uniform(r)*param.box[1];

        it = 0;
        //iteration over timesteps (corresponding to 1 experiment
        do
        {
        	cout<<i_iter<<" "<<it<<endl;
            i = 0;
            
            //loop over particles with positive energies
            do
            {
                x[0][i] = x[0][i] + sigma * rand_normal(0,1);
                x[1][i] = x[1][i] + sigma * rand_normal(0,1);

                BC(x[0][i],x[1][i],param);

                ix = NINT(x[0][i]/param.dx);  
                hist[it][ix] = hist[it][ix] + u[i];
                i++;

            }while(i<n_part/2);

            //loop over particles with negtive energies
            i = param.nsteps*param.n_part;
            do
            {
                x[0][i] = x[0][i] + sigma * rand_normal(0,1);
                x[1][i] = x[1][i] + sigma * rand_normal(0,1);

                BC(x[0][i],x[1][i],param);

                ix = NINT(x[0][i]/param.dx);
                hist[it][ix] = hist[it][ix] + u[i];
                i++;                

            }while(i<param.nsteps*param.n_part+n_part/2);

            n_part = n_part + 2*param.n_part;
            it++;

        }while(it<param.nsteps); 
        i_iter++;

    }

    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            hist[i][j] = hist[i][j] / (float)param.iter;

	std::ofstream fp1;
	fp1.open("/home/arjunbk/Documents/MASTER_PROJECT/Code/serial_energy_histogram.dat",std::fstream::out);
	fp1<<std::fixed;

	for (int i = 0; i < param.nsteps; ++i)
	 	{
	        for (int j = 0; j < (NINT(param.box[0]/param.dx)+1); ++j)
                fp1<<std::setprecision(5)<<float(i)<<"   "<<std::setprecision(5)<<float(j*param.dx - param.box[0]/2.0)<<"   "<<std::setprecision(5)<<float(hist[i][j])<<endl; 
	        fp1<<"\n";
	    }
	fp1.close();
	
	t2=clock();
	float diff = float(t2)-float(t1);
	cout<<"\nTotal run time: "<<(diff)/CLOCKS_PER_SEC<<" sec"<<endl;  
}



