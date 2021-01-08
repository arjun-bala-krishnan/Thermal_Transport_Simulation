#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
#include <cstdlib>
#include <iomanip> 
#include <cmath>
#include <time.h>
#include <mpi.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#define NINT(a) ((a) >= 0.0 ? (int)((a)+0.5) : (int)((a)-0.5))

const double PI  =3.141592653589793238463;

using std::cin;
using std::cout;
using std::endl;
using std::vector;
using std::ofstream;
using std::ifstream;
using std::fstream;
using std::setprecision;

struct material
{
	int num = 10;  //number of material 
	double pos[2][100]; //the centre of under 100 inclusions stored
	float diff = 1.0; //diffusion coefficient
    float radius = 0.5; //radius of inclusion
    float area_percent = 0.0; 
    float Rbd; //interfacial thermal resistance to the other material
    float rho; //density of the material
    float Cp; //specific heat capacity of material
    float v; // speed of sound in material
    float prob; //probability to the other material, matrix.prob = 4/(matrix.rho * matrix.Cp * matrix.v * matrix.Rbd), 
    			//incl.prob = matrix.Cf * pow((incl.radius + matrix.sd),3) - pow(incl.radius,3) * matrix.prob / (pow(incl.radius,3) - pow((incl.radius - incl.sd*matrix.sd),3));
    float Cf; //thermal equilibrium factor
    float sd; //standard deviation
};

struct syst
{
	long int nsteps = 1000; // random walk steps
	int n_part = 200;  // number of walkers 
	float time = 0.1; // time increment 
	float dx = 0.1; //scaling factor
	int iter = 4; //number of iteration (parallelised)
	float box[2] = {10.0,10.0}; //lengths of x and y
    float error = 0.001;
    bool incl = true; //true means inclusion present
    int dim = 2; //diamensions
    int attribute = 1; // attribute 0 means holllow and 1 means solid
    int collision = 0; // collision 0 means reflective and 1 means stochastic deflection
    float sigma; //param.sigma = sqrt(param.time*matrix.diff)
};

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

void PBC(double& x, double& y, syst param) // addittionaly bring PBC and RBC in an infinite loop
{
    if(x < 0.0)
        x = -x;

    else if(x > param.box[0])
        x = 2*param.box[0] - x;
  
    if(y  > param.box[1])
        y = y - param.box[1];

    else if(y < 0.0)
        y = y + param.box[1];
}

void check_inside_hollow(material incl, syst& param, double x1, double y1, double& x2, double& y2) 
{
    double x,y,x3,y3,p,q,q1,A,B,C,D,pdash,qdash;
    double x2dot,y2dot,x2dash,y2dash,distxx2dot,distx1x2dot,distx1x2dash,distpx3,distpx2dash;
    double m, m_line, m_refl, c, c_line, c_refl;
    double distxx1,distxx2,distxx3,distx1x2,distpx2,distpx1;
    double neigbour_list[2][100];
    double dx = 2*incl.radius; 
    double dy = 2*incl.radius;
    double weight;
	int ix,iy,jx,jy,kx,ky,flag,list_num;
	int ticket = 0;
	int check2_circle = 0;
	// int count = 0;


	if(x2 < 0.0 || y2< 0.0 || x2 > param.box[0] || y2 > param.box[1]) // Check boundary condition
		{
			PBC(x2,y2,param);
			ticket = 1; // only if the walker BC got changed ticket, it is allowed to enter BC (plus during one other special case)
		}
    do
    {
    	list_num = 0, x = 0.0, y = 0.0, x3 = 0.0, y3 = 0.0, x2dot = 0.0, y2dot = 0.0, x2dash = 0.0, y2dash = 0.0, weight = 0.0, p = 0.0, q = 0.0, q1 = 0.0;	 
    	distxx2dot = 0.0,distx1x2dot = 0.0,distx1x2dash = 0.0,distpx2 = 0.0; distpx3 = 0.0, distpx2dash = 0.0;
	    m = 0.0, m_line = 0.0, m_refl = 0.0, c = 0.0, c_line = 0.0, c_refl = 0.0;
	    A = 0.0, B = 0.0, C = 0.0, D = 0.0, pdash = 0.0, qdash = 0.0, flag = 0;

	    for(int i=0; i<incl.num; i++) // Create Neighbour List
	    {
			// int circle_num = 1;
			p = incl.pos[0][i];
			q = incl.pos[1][i];
			q1 = incl.pos[1][i];

			ix = int(x2/dx);
	   		iy = int(y2/dy);

			if(iy == param.box[1]-1) //iy is an integer number, the closest integer under top boundary will be (boundary - 1)
	   		{
	    		p = incl.pos[0][i];
	    		q1 = incl.pos[1][i]-dy;
	    		if(q < 0.0)
					q1 = q1 + param.box[1];
				iy = int(double(y2-dy)/double(dy));
	   		}

	   		if(iy == 0) //similarly iy is 0 for bottom boundary
	   		{
	    		p = incl.pos[0][i];
				q1 = incl.pos[1][i]+dy;
				if(q>param.box[1])
					q1 = q1 - param.box[1];
				iy = int(double(y2+dy)/double(dy));
	   		}

	   		
			jx = int(double(p)/double(dx));
			jy = int(double(q1)/double(dy));
			kx = jx-ix;
			ky = jy-iy;


			if( kx>=-2 && kx<=2 && ky>=-2 && ky<=2)
			{
				neigbour_list[0][list_num] = p;
				neigbour_list[1][list_num] = q;
				// cout<<"circle "<<circle_num<<": "<<p<<","<<q<<endl;
				// circle_num++;
				list_num++;
			}
			
	    } // End Neighbour list

	    // cout<<"got past neighbour list for x2, y2: "<<x2<<","<<y2<<endl;

	    if(list_num == 0) // there are no neighbours to the point which means we can exit the function
	    {
	    	flag = 0;
	    	break; //no neighbours in the list for the point
	    }

	    for(int i=0;i<list_num;i++)
	    {
			flag = 0;
			p = neigbour_list[0][i];
			q = neigbour_list[1][i];

			distpx2 = sqrt(pow((q-y2),2)+pow((p-x2),2));
			// cout<<"distance "<<distpx2<<endl;
			// cout<<"radius-0.01 "<<incl.radius-0.01<<endl;
			if (distpx2 < incl.radius-0.01) // see if the point 2 is inside the circle first
			{
				check2_circle = 1;
				if (q < incl.radius || q > (param.box[1]-incl.radius)) ticket = 1; // if circle crosses and point dont still make ticket to 1 so that it can get inside BC condition
			}
			else continue; // if point 2 is not inside circle move onto next inclusion

			// cout<<"check2_circle value: "<<check2_circle<<endl;
			
			// if point 2 crosses boundary or circle crosses or both together boundary conditions
			if ((y1 > (param.box[1]-2*incl.radius) && y2 < 2*incl.radius && ticket == 1) || ( y2 > (param.box[1]-2*incl.radius) && y1 < 2*incl.radius && ticket == 1) || (y1 > (param.box[1]-2*incl.radius) && q < 2*incl.radius) || ( q > (param.box[1]-2*incl.radius) && y1 < 2*incl.radius) )
			{
				// cout<<"got inside circle: "<<p<<","<<q<<" for BC in: "<<x2<<","<<y2<<endl;
				y2dash = y2; //in such cases that the point dont cross the boundary
				x2dash = x2;
				pdash = p; // in such cases that the circles dont cross the boundary
				qdash = q;
				
				if((y1 > (param.box[1]-2*incl.radius) && y2 < 2*incl.radius)) // point crosses upwards
				{
					// cout<<"point crosses upwards"<<endl;
					y2dash = param.box[1] + y2;
					x2dash = x2;
					m = (y2dash-y1)/(x2dash-x1);
					c = y1-m*x1;
					x2dot = (param.box[1] - c)/m; // dot means the point intersecting boundary (to calculate weight)
					y2dot = param.box[1];

					if (check2_circle == 1) // a circle that did not cross but point 2 was inside it after crossing
					{
						pdash = p;
						qdash = param.box[1] + q;
					}
				}

				if(( y2 > (param.box[1]-2*incl.radius) && y1 < 2*incl.radius)) // point crosses downwards
				{
					// cout<<"point crosses downwards"<<endl;
					y2dash = 0.0 - (param.box[1]-y2);
					x2dash = x2;
					m = (y2dash-y1)/(x2dash-x1);
					c = y1-m*x1;
					x2dot = (0.0 - c)/m; // dot means the point intersecting boundary (to calculate weight)
					y2dot = 0.0;

					if (check2_circle == 1) // a circle that did not cross but point 2 was inside it after crossing
					{
						pdash = p;
						qdash = 0.0 - (param.box[1] - q);
					}
				}

				if((y1 > (param.box[1]-2*incl.radius) && q < 2*incl.radius)) // circle crosses downwards (centre at bottom portion)
				{
					// cout<<"circle crosses downwards (centre at bottom portion)"<<endl;
					y2dash = param.box[1] + y2;
					x2dash = x2;
					m = (y2dash-y1)/(x2dash-x1);
					c = y1-m*x1;
					x2dot = (param.box[1] - c)/m; // dot means the point intersecting boundary (to calculate weight)
					y2dot = param.box[1];
					pdash = p;
					qdash = param.box[1] + q;
				}

				if(( q > (param.box[1]-2*incl.radius) && y1 < 2*incl.radius)) // circle crosses upwards (centre on top portion)
				{
					// cout<<"circle crosses upwards (centre on top portion)"<<endl;
					y2dash = 0.0 - (param.box[1]-y2);
					x2dash = x2;
					m = (y2dash-y1)/(x2dash-x1);
					c = y1-m*x1;
					x2dot = (0.0 - c)/m; // dot means the point intersecting boundary (to calculate weight)
					y2dot = 0.0;
					pdash = p;
					qdash = 0.0 - (param.box[1] - q);
				}

				if ((qdash==q && y2dash!=y2 || qdash!=q && y2dash==y2) && check2_circle == 1) 
				{
					if (qdash==q && y2dash!=y2) // point crosses circle dont, but point is inside circle 
					{
						if (q > 0.5*param.box[1])
						{
							pdash = p;
							qdash = 0.0 - (param.box[1] - q);
						} 

						if (q < 0.5*param.box[1])
						{
							pdash = p;
							qdash = param.box[1] + q;
						}
					}

					if (qdash!=q && y2dash==y2) //circle crosses but point dont, but point is inside circle
					{
						if (y2 > 0.5*param.box[1])
						{
							y2dash = param.box[1] + y2;
							x2dash = x2;
						} 

						if (y2 < 0.5*param.box[1])
						{
							y2dash = 0.0 - (param.box[1]-y2);
							x2dash = x2;
						}
					}
				}

				distpx2dash = sqrt(pow((pdash-x2dash),2)+pow((qdash-y2dash),2)); 
				if(distpx2dash>incl.radius) continue; //check if the virtual point 2 is outside the circle, check and continue only if the 2nd point is in the circle    
	
				m = (y2dash-y1)/(x2dash-x1);
				c = y1 - m*x1;				
				A = pow(m,2)+1.0;
				B = 2*(m*c - m*qdash  - pdash);
				C = (pow(qdash,2) - pow(incl.radius,2) + pow(pdash,2) - 2*c*qdash  + pow(c,2));
				D = (pow(B,2)-4*A*C);

				if(D < 0.0) continue; // circle does not intersect line, so continue to next point 

				x = ( -1*B + sqrt(D) )/ ( 2*A ); 
				y = m*x + c;

				distxx1 = sqrt(pow((y-y1),2)+pow((x-x1),2));
				distx1x2dash = sqrt(pow((y1-y2dash),2)+pow((x1-x2dash),2));
			
				if(distx1x2dash<distxx1) // checking if intercept is between start and end points
				{
					x = ( -1*B - sqrt(pow(B,2)-4*A*C) )/ ( 2*A ); 
					y = m*x + c;
					distxx1 = sqrt(pow((y-y1),2)+pow((x-x1),2));
				}  

				if((y1 > (param.box[1]-incl.radius) && y2 < incl.radius)) // calculating weight if point crossed downwards 
				{
					distxx2dot = sqrt(pow((y-y2dot),2)+pow((x-x2dot),2));
					weight = distxx2dot + sqrt(pow((y2-0.0),2)+pow((x2-x2dot),2));
				}
					
				if(( y2 > (param.box[1]-incl.radius) && y1 < incl.radius)) // calculating weight if point crossed upwards
				{
					distxx2dot = sqrt(pow((y-y2dot),2)+pow((x-x2dot),2));
					weight = distxx2dot + sqrt(pow((param.box[1]-y2),2)+pow((x2-x2dot),2));
				}

				m_line = -1.0/((qdash-y)/(pdash-x));  //contructing tangent line
				c_line = y - m_line*x;             
				m_refl = (2*m_line + m*pow(m_line,2) - m) / (2*m_line*m - pow(m_line,2) +1.0); //constructing reflection line
				c_refl = y - m_refl*x;

			} // so now we have x, y ,m ,c , tangent, reflection and weight

			else // start of normal cases there is no dot or dash values
			{
				// cout<<"got inside circle: "<<p<<","<<q<<" for nrmal condition in: "<<x2<<","<<y2<<endl;
				m = (y2-y1)/(x2-x1);
				c = y1 - m*x1;
				A = pow(m,2)+1.0;
				B = 2*(m*c - m*q - p);
				C = (pow(q,2) - pow(incl.radius,2) + pow(p,2) - 2*c*q + pow(c,2));
				D = (pow(B,2)-4*A*C);
	
				if(D < 0.0) continue; //circle does not intersect line, so continue to next point

				x = ( -1*B + sqrt(D) )/ ( 2*A ); 
				y = m*x + c;
				distxx1 = sqrt(pow((y-y1),2)+pow((x-x1),2));
				distxx2 = sqrt(pow((y-y2),2)+pow((x-x2),2));
				distx1x2 = sqrt(pow((y2-y1),2)+pow((x2-x1),2)); 

				if(distx1x2<distxx1) // checking if intercept is between start and end points
				{
					x = ( -1*B - sqrt(pow(B,2)-4*A*C) )/ ( 2*A ); 
					y = m*x + c;
					distxx1 = sqrt(pow((y-y1),2)+pow((x-x1),2));
					distxx2 = sqrt(pow((y-y2),2)+pow((x-x2),2));  
				}  

				distpx2 = sqrt(pow((p-x2),2)+pow((q-y2),2)); 
				if(distpx2>incl.radius) continue;
				weight = distxx2; // assigning weight
				m_line = -1.0/((q-y)/(p-x));  //contructing tangent line
				c_line = y - m_line*x;             
				m_refl = (2*m_line + m*pow(m_line,2) - m) / (2*m_line*m - pow(m_line,2) +1.0); //constructing reflection line
				c_refl = y - m_refl*x;
				
			} // so now we have x, y ,m ,c , tangent, reflection and weight

	        if(param.collision == 1) // stochastic dispersion
	        {
	        	float rand_fract = 2.0*rand()/RAND_MAX - 1;
	        	m_refl = 100.0 * rand_fract; //random value of slope with maximum 100 and minimum -100 slope
	        	c_refl = y - m_refl*x;
	        }

	        x3 = x + sqrt(pow(weight,2)/(1.0+pow(m_refl,2))); // constructing reflection point
	        y3 = m_refl*(x3-x) + y;

	        distpx3 = pow((x3-p),2) + pow((y3-q),2); // to compare square of distance in normal cases
	        if((y1 > (param.box[1]-incl.radius) && q < incl.radius) || ( q > (param.box[1]-incl.radius) && y1 < incl.radius)) 
				distpx3 = pow((x3-pdash),2) + pow((y3-qdash),2); // square of distance circle cross case
			
			if ( distpx3 < pow(incl.radius,2) ) // if weighted x3 is inside inclusion then we have to throw it outward
			{
				x3 = x - sqrt(pow(weight,2)/(1.0+pow(m_refl,2)));
				y3 = m_refl*(x3-x) + y;
			}
				
			if(x3 < 0.0 || y3< 0.0 || x3 > param.box[0] || y3 > param.box[1]) // Check boundary condition for reflected point
				PBC(x3,y3,param);
			
			x2=x3;
			y2=y3;
			flag = 1;
			break;

		} // End of loop inside neghbouring list

    }while (flag == 1); //check the new point from beginning till end, starting from neighbour list creation

} // End of check_inside_hollow fuction

void check_inside_solid(material incl, material matrix, syst& param, double x1, double y1, double& x2, double& y2)
{
	double x,y,x3,y3,p,q,A,B,C,D,m,c,weight,x2x1;
    double distxx1,distxx2,distxx3,distx1x2,distpx2,distpx1,distpx3,distx3x2;
    double dx = 2*incl.radius; 
    double dy = 2*incl.radius;
    double neigbour_list[2][100];
    double diff_dist = 0.0, distx1x3 = 0.0;
    double diff_fract = sqrt(incl.diff/matrix.diff);
	int ticket = 0, ticket2 = 0; //start ticket = 0 
	int check1_circle = 0, check2_circle = 0;
	int indx = 0;

	// cout<<"x1 and y1 during function start: "<<x1<<","<<y1<<endl;

	int ix,iy,jx,jy,kx,ky,flag,list_num = 0;
	ix = 0, iy = 0, jx = 0, jy = 0, kx = 0, ky = 0;

	for(int i=0; i<incl.num; i++) // Create Neighbour List for point 1
    {
		p = incl.pos[0][i];
		q = incl.pos[1][i];
		ix = int(double(x1)/double(dx));
   		iy = int(double(y1)/double(dy));

		if(iy == (param.box[1]-1) && q <= dy) q = param.box[1]+q; // iy is an integer number, the closest integer under top boundary will be (boundary - 1)
   		if(iy == 0 && q > (param.box[1]-dy)) q = 0.0 - (param.box[1]-q);// similarly iy is 0 for bottom boundary

		jx = int(double(p)/double(dx));
		jy = int(double(q)/double(dy));
		kx = jx-ix;
		ky = jy-iy;



		if( (kx<=1 && kx>=-1) && (ky<=1 && ky>=-1) )
		{
			neigbour_list[0][list_num] = p;
			neigbour_list[1][list_num] = q;
			list_num++;
		}

	} // End Neighbour list

	for (int i=0;i<list_num;i++)
	{
		p = neigbour_list[0][i];
		q = neigbour_list[1][i];

		if(y1>(param.box[1]-dy) && q<dy) q = param.box[1]+q;
		if(y1<dy && q>(param.box[1]-dy)) q = 0.0 - (param.box[1]-q);

		distpx1 = sqrt(pow((q-y1),2)+pow((p-x1),2));

		if (distpx1 < incl.radius-0.01) // check if the 1st point is inside the inclusion
			{
				check1_circle = 1;
				indx = i;
				break;
			}
	}

	if(check1_circle == 1)
	{
		p = neigbour_list[0][indx];
		q = neigbour_list[1][indx]; 
		x2 = x2;
		y2 = y2; 
		x = 0.0;
		y = 0.0;
		check2_circle = 0;

		m = (y2-y1)/(x2-x1);
		c = y1 - m*x1;
		A = pow(m,2)+1.0;
		B = 2*(m*c - m*q - p);
		C = (pow(q,2) - pow(incl.radius,2) + pow(p,2) - 2*c*q + pow(c,2));
		D = (pow(B,2)-4*A*C);;

		if(D < 0.0) D = 0.0; // error as the circle definitly meets the line as point 1 is inside the circle

		x = ( -1*B + sqrt(D) )/ ( 2*A ); 
		y = m*x + c;

		distxx2 = sqrt(pow((y2-y),2)+pow((x2-x),2)); 
		distx1x2 = sqrt(pow((y2-y1),2)+pow((x2-x1),2)); 
		distxx1 = sqrt(pow((y1-y),2)+pow((x1-x),2)); 

		distpx2 = sqrt(pow((y2-q),2)+pow((x2-p),2)); 
		if (distpx2 < incl.radius ) check2_circle = 1; // This means that point 2 is inside the circle

		if(((distxx1 + distxx2 - distx1x2 > 0.01) && check2_circle == 0) || ((distx1x2 + distxx2 - distxx1 > 0.01) && check2_circle == 1)) // checking if intercept is between start and end points or if end point is inside inclusion is intercept in the same direction
		{
			x = ( -1*B - sqrt(pow(B,2)-4*A*C) )/ ( 2*A ); 
			y = m*x + c;
			distxx1 = sqrt(pow((y-y1),2)+pow((x-x1),2));
			distxx2 = sqrt(pow((y2-y),2)+pow((x2-x),2));
		}  

		distx1x3 = 0.0;

		if (distxx1 + distxx2 - distx1x2 < 0.1) // condition when point 2 is outside the circle
		{
			weight = distxx1 * diff_fract;
			diff_dist = distxx1 - weight;
			distx1x3 = distx1x2 - diff_dist;
		}

		if ((distx1x2 + distxx2 - distxx1 < 0.01) || (distx1x2 - distxx1 < 0.01)) // condition when point 2 is inside or on the circle
			distx1x3 = distx1x2 * diff_fract;
		
		x3 = x1 + sqrt(pow(distx1x3,2)/(1.0+pow(m,2))); // constructing new point 2 as 3
    	y3 = m*x3 + c;
    	distx3x2 = sqrt(pow((y2-y3),2)+pow((x2-x3),2));    

    	if (distx1x3 + distx3x2 - distx1x2 >0.01) // in both case 1 and case 2, this condition should not happen, so check and change if required
    	{
    		x3 = x1 - sqrt(pow(distx1x3,2)/(1.0+pow(m,2))); // constructing new point 2
    		y3 = m*x3 + c;
    	}

    	distpx3 = sqrt(pow((q-y3),2)+pow((p-x3),2));
 
	    if (x3 < 0.0 || y3< 0.0 || x3 > param.box[0] || y3 > param.box[1]) 
	    {
	    	PBC(x3,y3,param); // Check boundary condition for new point 3
	    	ticket = 1;
	    }

		x2 = x3;
		y2 = y3;	

		if (distpx3 < incl.radius + 0.01) 
		{	
			PBC(x2,y2,param);			
			return; // the point 2 is still inside the inclusion stil has to be returned (because no inclusions overlap allowed)	
		}
    }
    
	list_num = 0, x = 0.0, y = 0.0, x3 = 0.0, y3 = 0.0, weight = 0.0, p = 0.0, q = 0.0;
	ix = 0, iy = 0, jx = 0, jy = 0, kx = 0, ky = 0;	 
	distx1x2 = 0.0,distpx2 = 0.0; distpx3 = 0.0, distpx2 = 0.0;
    m = 0.0, c = 0.0;
    A = 0.0, B = 0.0, C = 0.0, D = 0.0, flag = 0;

    for(int i=0; i<incl.num; i++) // Create Neighbour List for point 2
    {
		p = incl.pos[0][i];
		q = incl.pos[1][i];
		ix = int(double(x2)/double(dx));
   		iy = int(double(y2)/double(dy));

		if(iy >= (param.box[1]-1) && q <= dy) q = param.box[1]+q; // iy is an integer number, the closest integer under top boundary will be (boundary - 1)
   		if(iy <= 0 && q > (param.box[1]-dy)) q = 0.0 - (param.box[1]-q);// similarly iy is 0 for bottom boundary

		jx = int(double(p)/double(dx));
		jy = int(double(q)/double(dy));
		kx = jx-ix;
		ky = jy-iy;

		if( (kx<=1 && kx>=-1) && (ky<=1 && ky>=1) )
		{
			neigbour_list[0][list_num] = p;
			neigbour_list[1][list_num] = q;
			list_num++;
		}
	} // End Neighbour list


    if(list_num == 0) // for cases where point 2 is not inside any inclusion
    {	
		if (x2 < 0.0 || y2< 0.0 || x2 > param.box[0] || y2 > param.box[1]) PBC(x2,y2,param); // Check boundary condition for point 2
		return;// There are no neighbours to the point so the point 2 can be returned after BC
    }


    for(int i=0;i<list_num;i++)
    {
		p = neigbour_list[0][i];
		q = neigbour_list[1][i];
		check2_circle = 0;

		if( y2>(param.box[1]-dy) && q<dy ) q = (param.box[1]+q);
		if( y2<dy && q>(param.box[1]-dy) ) q = 0.0 - (param.box[1]-q);

		distpx2 = sqrt(pow((q-y2),2)+pow((p-x2),2));
		if (distpx2 < incl.radius-0.01) // check if the 1st point is inside the inclusion
			{
				check2_circle = 1;
				indx = i;
			}
		
		if(check2_circle == 0) continue; // if circle is not inside the inclusion

		x2x1 = (x2-x1);
		if (fabs(x2x1) < fabs(0.001)) x2x1 = 0.001; //minor error is rejected in case of vertical lines

		m = (y2-y1)/(x2x1);
		c = y1 - m*x1;
		A = pow(m,2)+1.0;
		B = 2*(m*c - m*q - p);
		C = (pow(q,2) - pow(incl.radius,2) + pow(p,2) - 2*c*q + pow(c,2));
		D = (pow(B,2)-4*A*C);

		if(D < 0.0)  // D must be always positive
		{
			cout<<"D error in non BC condition"<<endl;
			if(x2 < 0.0 || y2< 0.0 || x2 > param.box[0] || y2 > param.box[1]) PBC(x2,y2,param);
			return;
		}

		x = ( -1*B + sqrt(D) )/ ( 2*A ); 
		y = m*x + c;
		distxx1 = sqrt(pow((y-y1),2)+pow((x-x1),2));
		distxx2 = sqrt(pow((y-y2),2)+pow((x-x2),2));
		distx1x2 = sqrt(pow((y2-y1),2)+pow((x2-x1),2)); 

		if(distx1x2< (distxx1)) // checking if intercept is between start and end points
		{
			x = ( -1*B - sqrt(pow(B,2)-4*A*C) )/ ( 2*A ); 
			y = m*x + c;
			distxx1 = sqrt(pow((y-y1),2)+pow((x-x1),2));
			distxx2 = sqrt(pow((y-y2),2)+pow((x-x2),2));  
		} 

		weight = distxx2; // assigning weight		
		weight = weight * diff_fract;
		x3 = x + sqrt(pow(weight,2)/(1.0+pow(m,2))); // constructing transmission point
	    y3 = m*(x3-x) + y; 

	    distpx3 = sqrt(pow((x3-p),2) + pow((y3-q),2)); 
		
		if ( distpx3 > incl.radius ) // to make sure transmitted x3 and y3 is inside inclusion
		{
			x3 = x - sqrt(pow(weight,2)/(1.0+pow(m,2))); // constructing transmission point
	    	y3 = m*(x3-x) + y;
	    	distpx3 = sqrt(pow((x3-p),2) + pow((y3-q),2)); 
	    }	 

		if(x3 < 0.0 || y3< 0.0 || x3 > param.box[0] || y3 > param.box[1]) PBC(x3,y3,param); // Check boundary condition for transmitted point
		x2=x3;
		y2=y3;
		break; // no longer need to continue as it is supposed to stay inside the current inclusion	

	} // End of loop inside neghbouring list

	if(x2 < 0.0 || y2< 0.0 || x2 > param.box[0] || y2 > param.box[1]) PBC(x2,y2,param); // Check boundary condition for transmitted point

} // End of fuction

int main(int argc, char** argv) 
{
    syst param;
    srand(time(NULL));
    int i,k,ix,iy,it_steps,n_part;
    const int rows = param.nsteps;
    const int cols = NINT(param.box[0]/param.dx)+1;
    float hist[rows][cols];
    float hist_avg[rows][cols];
    double x1,y1;
    material matrix; 
    material incl;
    matrix.diff = 0.712; 
    incl.diff = 0.135; 
    matrix.num = 1;
    param.sigma = sqrt(param.time*matrix.diff);    
    matrix.sd = 1.0;
    
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
        {
            hist[i][j] = 0.0;
            hist_avg[i][j] = 0.0;
        }

    double *u = (double*) malloc(sizeof(double)*(param.n_part*2*param.nsteps)); 
    
    for(int i=0;i<(param.nsteps*param.n_part);i++)
        u[i]=1.0;

    for(int i=(param.nsteps*param.n_part);i<param.nsteps*2*param.n_part;i++)
        u[i]=-1.0;

    MPI_Init( &argc,&argv);
    double tbeg = MPI_Wtime();  
    int rank,size,steps,indx,start,stop;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int count = param.iter/size;
	int reminder = param.iter % size;

	if(rank == 0)
	{
		gsl_rng * r;
        const gsl_rng_type * T;
        T = gsl_rng_gfsr4;
        r = gsl_rng_alloc (T);
        gsl_rng_set (r, rand());		
		double x_y,ran;
        int flag = 1;

        if(param.incl == true) // making random inclusions
        {
			if(incl.area_percent == 0)
        		incl.area_percent = (PI*pow(incl.radius,2) * incl.num * 100.0) / (param.box[0]*param.box[1]);
			else
				incl.num = int(incl.area_percent * (param.box[0] * param.box[1]) / (PI * pow(incl.radius,2) * 100.0));

			cout<<"No: of inclusions are: "<<incl.num<<"\nPercentage area: "<<incl.area_percent<<endl;
			double dist;
			double x0,y0,x0dash,y0dash;

			for(int i=0; i<incl.num; i++) 
				{	
					flag = 0;
					while(1)
					{						

						for (int j=0; j<2; j++) //createx and y for the ith inclusion
						{	
							ran = gsl_rng_uniform(r);
							incl.pos[j][i]= ran*param.box[j]; 
						}

						x0 = incl.pos[0][i]; // save them as x0, y0
						y0 = incl.pos[1][i];

						if( (x0 < 4*incl.radius) || (x0 > (param.box[0]-4*incl.radius)) ) flag = 1; //check if x0 values satisfy boundary condition

						if(flag==1) 
							{
								i--; 
								break;
							}

						for (int j = 0; j < i ; j++) // check distance is below diameter

						{
							dist = sqrt(pow((incl.pos[0][j]-x0),2) + pow((incl.pos[1][j]-y0),2));

							if((incl.pos[1][j]>(param.box[1]-incl.radius) && y0<incl.radius)) 
							{
								x0dash = x0;
								y0dash = param.box[1] + y0;
								dist = sqrt(pow((incl.pos[0][j]-x0dash),2) + pow((incl.pos[1][j]-y0dash),2));
							}

							if((y0>(param.box[1]-incl.radius) && incl.pos[1][j]<incl.radius))
							{
								x0dash = x0;
								y0dash = 0.0 - (param.box[1] - y0);
								dist = sqrt(pow((incl.pos[0][j]-x0dash),2) + pow((incl.pos[1][j]-y0dash),2));
							}

							if(dist < 2*incl.radius) 
								{
									flag = 1;
									break;
								}			 
						}

						if(flag==1) 
							{
								i--; 
								break;
							}

			            if( flag == 0 ) break;            // move onto next cordinate

					}
				}	
		
			std::ofstream fp1;
			fp1.open("/home/arjunbk/Documents/MASTER_PROJECT/Code/inclusions.dat",std::fstream::out);        
	    	fp1<<std::fixed;
			cout<<"Inclusion cordinates are: "<<endl;
	        for(int i=0; i<incl.num; i++) 
	        {
	        	cout<<"("<<i+1<<") "<<incl.pos[0][i]<<", "<<incl.pos[1][i]<<endl;
	        	fp1<<i+1<<" "<<std::setprecision(5)<<float(incl.pos[0][i])<<"  "<<std::setprecision(5)<<float(incl.pos[1][i])<<" "<<std::setprecision(5)<<float(incl.radius)<<endl;
        		if(incl.pos[1][i] < incl.radius) 
				{
					x0dash = incl.pos[0][i];
					y0dash = param.box[1] + incl.pos[1][i];
					fp1<<i+1<<"' "<<std::setprecision(5)<<float(x0dash)<<"  "<<std::setprecision(5)<<float(y0dash)<<" "<<std::setprecision(5)<<float(incl.radius)<<endl;
				}

				if(incl.pos[1][i] > (param.box[1]-incl.radius)) 
				{
					x0dash = incl.pos[0][i];
					y0dash = 0.0 - (param.box[1] - incl.pos[1][i]);
					fp1<<i+1<<"' "<<std::setprecision(5)<<float(x0dash)<<"  "<<std::setprecision(5)<<float(y0dash)<<" "<<std::setprecision(5)<<float(incl.radius)<<endl;
				}

	        }
		    cout<<endl;	
		    fp1.close();	
		}  	
	}

	MPI_Bcast(incl.pos, 100*2 , MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	if (rank < reminder) // The first 'remainder' ranks get 'count + 1' tasks each
	{	    
	    start = rank * (count + 1);
	    stop = start + count;
	} 
	else // The remaining 'size - remainder' ranks get 'count' task each
	{
	    start = rank * count + reminder;
	    stop = start + (count - 1);
	}

	for (int itr = start; itr <= stop; ++itr)  
	{	
		gsl_rng * r;
        const gsl_rng_type * T;
        T = gsl_rng_gfsr4;
        r = gsl_rng_alloc (T);
        gsl_rng_set (r, (itr+1)*rand());        
		
	    double **x = (double**)malloc(sizeof(double*) * 2);
	    for (int i = 0; i < 2; i++) 
	        x[i] = (double *)malloc(sizeof(double) * (param.n_part*2*param.nsteps));

	    for(int i=0;i<(param.nsteps*param.n_part);i++) 
	        x[0][i] = 0.0;

	    for(int i=(param.nsteps*param.n_part);i<2*param.nsteps*param.n_part;i++) 
	        x[0][i] = param.box[0];

	    for(int i=0;i<2*param.nsteps*param.n_part;i++) 
	        x[1][i] = gsl_rng_uniform(r)*param.box[1];
 		
 		n_part = 2*param.n_part;; //the number of walkers added at both ends
        it_steps = 0;     

        do
        {
	        i = 0;
	        cout<<it_steps<<" "<<itr<<" "<<endl;
	        x1 = x[0][i];
	       	y1 = x[1][i];
            do
            {	           	
            	if(i>0)
                {
                	x1 = x[0][i];
	       			y1 = x[1][i];
                }

            	x[0][i] = x[0][i] + param.sigma * rand_normal(0,1);
                x[1][i] = x[1][i] + param.sigma * rand_normal(0,1);      

     			if(param.incl==true)
     			{
     				if (param.attribute == 0)
						check_inside_hollow(incl,param,x1,y1,x[0][i],x[1][i]);
					if (param.attribute == 1)
						check_inside_solid(incl,matrix,param,x1,y1,x[0][i],x[1][i]);
     			}

				else    
					PBC(x[0][i],x[1][i],param);  
				
				ix = NINT(x[0][i]/param.dx); 
                hist[it_steps][ix] = hist[it_steps][ix] + u[i];
                x1 = x[0][i];
	       		y1 = x[1][i];
                i++;  

            }while(i<n_part/2);

            i = param.nsteps*param.n_part;
            x1 = x[0][i];
           	y1 = x[1][i];           

            do
            {                
                if(i>param.nsteps*param.n_part)
                {
                	x1 = x[0][i];
	       			y1 = x[1][i];
                }

                x[0][i] = x[0][i] + param.sigma * rand_normal(0,1);
                x[1][i] = x[1][i] + param.sigma * rand_normal(0,1);

                if(param.incl==true)
                {
					if (param.attribute == 0)
						check_inside_hollow(incl,param,x1,y1,x[0][i],x[1][i]); 
					if (param.attribute == 1)
						check_inside_solid(incl,matrix,param,x1,y1,x[0][i],x[1][i]); 
                }
				
				else    
					PBC(x[0][i],x[1][i],param);

				ix = NINT(x[0][i]/param.dx);
                hist[it_steps][ix] = hist[it_steps][ix] + u[i];
                x1 = x[0][i];
                y1 = x[1][i];
                i++; 

            }while(i<param.nsteps*param.n_part + n_part/2);
         
            n_part = n_part + 2*param.n_part;
            it_steps++;

        }while(it_steps<param.nsteps);
	}

	MPI_Barrier( MPI_COMM_WORLD );
	MPI_Reduce( &hist, &hist_avg, rows*cols, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD );

	for (int i = 0; i < param.nsteps; ++i)
    	for (int j = 0; j < (NINT(param.box[0]/param.dx)+1); ++j)
        	hist_avg[i][j] = hist_avg[i][j]/param.iter;
                
    if ( rank == 0 ) 
    {
        std::ofstream fp1,fp2;
	    fp1.open("/home/arjunbk/Documents/MASTER_PROJECT/Code/parallel_energy_histogram.dat",std::fstream::out);
        fp2.open("/home/arjunbk/Documents/MASTER_PROJECT/Code/convergence.dat",std::fstream::out);
        
	    fp1<<std::fixed;
        fp2<<std::fixed;

        float a,b;
        float *y = (float*) malloc(sizeof(float)*(param.nsteps)); 
        for(int i=0;i<param.nsteps;i++)
        y[i]=0.0;
        float *dy = (float*) malloc(sizeof(float)*(param.nsteps)); 
        for(int i=0;i<param.nsteps;i++)
        dy[i]=0.0;

		for (int i = 0; i < param.nsteps; ++i)
		 	{
                for (int j = 0; j < (NINT(param.box[0]/param.dx)+1); ++j)
                    {
                        fp1<<std::setprecision(5)<<float(i)<<"   "<<std::setprecision(5)<<float(j*param.dx - param.box[0]/2.0)
                        											<<"   "<<std::setprecision(5)<<float(hist_avg[i][j])<<endl; 

                        a = 0.0;
                        b = 0.0;
                        if(j%param.nsteps==0)
                        {
                            for (int k = j; k < j+param.nsteps; ++k)
                                y[i] = y[i]+ (hist[i][j] - hist[i][j+99]);

                            y[i] = y[i]/float(param.nsteps);
                        }
                    }
                fp1<<"\n";  
            }
            int flag = 0;
            for(int i=0;i<param.nsteps-1;i++)
            {
                dy[i] = y[i+1]-y[i];
                dy[i] = dy[i]/i;
                fp2<<std::setprecision(5)<<float(i)<<"      "<<std::setprecision(5)<<dy[i]<<endl;
                if (dy[i]>-param.error && dy[i]<param.error && flag==0)
                    {
                        flag = 1;
                        cout<<"convergence achieved at step: "<<i<<" with value = "<<dy[i]<<endl;
                    }
            }
            if(flag==0)
                cout<<"convergence failed to achieve"<<endl;

	    fp1.close();
	    fp2.close();      
    }

    MPI_Barrier( MPI_COMM_WORLD );
    double elapsedTime = MPI_Wtime() - tbeg;
    double totalTime;
    MPI_Reduce( &elapsedTime, &totalTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD ); 

    if(rank == 0)
        cout<<"\nTotal time spent in seconds is: "<<totalTime<<endl;  

  	MPI_Finalize();
}
