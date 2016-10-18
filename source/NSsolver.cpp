#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <vector>
#include "init.h"
#include "gridgen.h"
#include "metric.h"
#include "RHS.h"
#include "BC.h"


using namespace std;

int main()
{
	
	int N= 61;                          /*Grid size`*/
	double st=1; 			  /*Stretching ratio, change this according to need*/
	double Re=100; 		/*Reynolds number*/
	double eps=0.01;   /*dissipation factor*/
	double ep=0.01;		/*residual smoothing factor*/    
	double L=30; /*Length of the domain*/
	double t=0;
	double t1=0.1;		/*Time duration, in seconds, for the study. Put time according to need or put a very large number for steady state problems*/
	double dt=0.005;	/*Time step*/
	double xmax=2*M_PI;	/* range of the grid in x-direction, change this according to the problem*/
	double ymax=2*M_PI; /* range of the grid in y-direction, change this according to the problem*/
	vector<vector<double> > x (N,vector<double>(N, 0));
	vector<vector<double> > y (N,vector<double>(N, 0));
	vector<vector<double> > xvel (N,vector<double>(N, 0));	
	vector<vector<double> > yvel (N,vector<double>(N, 0));
	vector<vector<double> > xvel1 (N,vector<double>(N, 0));
	vector<vector<double> > yvel1 (N,vector<double>(N, 0));
	vector<vector<double> > Press (N,vector<double>(N, 0));
	
	
	gridgen(N, st, xmax, ymax, x, y);  // Grid generation
	
	
	init( N, Re, t, t1, x, y, xvel, xvel1, yvel, yvel1, Press); // Initial conditions
	
	vector<vector<double> > u_new (N,vector<double>(N, 0));
	vector<vector<double> > u_new1 (N,vector<double>(N, 0));
	vector<vector<double> > u_new2 (N,vector<double>(N, 0));
	vector<vector<double> > u_new3 (N,vector<double>(N, 0));
	vector<vector<double> > u_k (N,vector<double>(N, 0));
	vector<vector<double> > u_k1 (N,vector<double>(N, 0));
	vector<vector<double> > u_k2 (N,vector<double>(N, 0));
	vector<vector<double> > u_k3 (N,vector<double>(N, 0));
	vector<vector<double> > u_n (N,vector<double>(N, 0));
	vector<vector<double> > u_n1 (N,vector<double>(N, 0));
	vector<vector<double> > u_n2 (N,vector<double>(N, 0));
	vector<vector<double> > u_n3 (N,vector<double>(N, 0));
	vector<vector<double> > u_old (N,vector<double>(N, 0));
	vector<vector<double> > u_old1 (N,vector<double>(N, 0));
	vector<vector<double> > u_old2 (N,vector<double>(N, 0));
	vector<vector<double> > u_old3 (N,vector<double>(N, 0));

	vector<vector<double> > v_new (N,vector<double>(N, 0));
	vector<vector<double> > v_new1 (N,vector<double>(N, 0));
	vector<vector<double> > v_new2 (N,vector<double>(N, 0));
	vector<vector<double> > v_new3 (N,vector<double>(N, 0));
	vector<vector<double> > v_k (N,vector<double>(N, 0));
	vector<vector<double> > v_k1 (N,vector<double>(N, 0));
	vector<vector<double> > v_k2 (N,vector<double>(N, 0));
	vector<vector<double> > v_k3 (N,vector<double>(N, 0));
	vector<vector<double> > v_n (N,vector<double>(N, 0));
	vector<vector<double> > v_n1 (N,vector<double>(N, 0));
	vector<vector<double> > v_n2 (N,vector<double>(N, 0));
	vector<vector<double> > v_n3 (N,vector<double>(N, 0));
	vector<vector<double> > v_old (N,vector<double>(N, 0));
	vector<vector<double> > v_old1 (N,vector<double>(N, 0));
	vector<vector<double> > v_old2 (N,vector<double>(N, 0));
	vector<vector<double> > v_old3 (N,vector<double>(N, 0));

	vector<vector<double> > p_new (N,vector<double>(N, 0));
	vector<vector<double> > p_new1 (N,vector<double>(N, 0));
	vector<vector<double> > p_new2 (N,vector<double>(N, 0));
	vector<vector<double> > p_new3 (N,vector<double>(N, 0));
	vector<vector<double> > p_k (N,vector<double>(N, 0));
	vector<vector<double> > p_k1 (N,vector<double>(N, 0));
	vector<vector<double> > p_k2 (N,vector<double>(N, 0));
	vector<vector<double> > p_k3 (N,vector<double>(N, 0));
	vector<vector<double> > p_n (N,vector<double>(N, 0));
	vector<vector<double> > p_n1 (N,vector<double>(N, 0));
	vector<vector<double> > p_n2 (N,vector<double>(N, 0));
	vector<vector<double> > p_n3 (N,vector<double>(N, 0));
	vector<vector<double> > p_old (N,vector<double>(N, 0));
	vector<vector<double> > p_old1 (N,vector<double>(N, 0));
	vector<vector<double> > p_old2 (N,vector<double>(N, 0));
	vector<vector<double> > p_old3 (N,vector<double>(N, 0));

	vector<vector<double> > dtau (N,vector<double>(N, 0));
	vector<vector<double> > JC (N,vector<double>(N, 0));
	vector<vector<double> > ex (N,vector<double>(N, 0));
	vector<vector<double> > ey (N,vector<double>(N, 0));
	vector<vector<double> > zx (N,vector<double>(N, 0));
	vector<vector<double> > zy (N,vector<double>(N, 0));
	vector<vector<double> > rus (N,vector<double>(N, 0));
	vector<vector<double> > rvs (N,vector<double>(N, 0));
	vector<vector<double> > rcs (N,vector<double>(N, 0));
	vector<vector<double> > rho1 (N,vector<double>(N, 0));
	vector<vector<double> > rho2(N,vector<double>(N, 0));

	vector<vector<double> > RHSu11(N,vector<double>(N, 0));
	vector<vector<double> > RHSu22(N,vector<double>(N, 0));
	vector<vector<double> > RHSu33(N,vector<double>(N, 0));
	vector<vector<double> > RHSv11(N,vector<double>(N, 0));
	vector<vector<double> > RHSv22(N,vector<double>(N, 0));
	vector<vector<double> > RHSv33(N,vector<double>(N, 0));
	vector<vector<double> > RHSu(N,vector<double>(N, 0));
	vector<vector<double> > RHSv(N,vector<double>(N, 0));


	metric( N, x, y,zx,zy,ex,ey,JC);   //Metric Calculation

	for (int i=0; i<N; i++)  // Initializing all the velocity variables and dtau
	{
		for (int j=0; j<N; j++)
		{
		
			u_new1[i][j]=xvel[i][j];
			u_new2[i][j]=xvel[i][j];
			u_new3[i][j]=xvel[i][j];
			u_new[i][j]=xvel[i][j];
			u_k1[i][j]=xvel[i][j];
			u_k2[i][j]=xvel[i][j];
			u_k3[i][j]=xvel[i][j];
			u_k[i][j]=xvel[i][j];
			u_n1[i][j]=xvel[i][j];
			u_n2[i][j]=xvel[i][j];
			u_n3[i][j]=xvel[i][j];
			u_n[i][j]=xvel[i][j];
			u_old1[i][j]=xvel[i][j];
			u_old2[i][j]=xvel[i][j];
			u_old3[i][j]=xvel[i][j];
			u_old[i][j]=xvel[i][j];
			v_new1[i][j]=yvel[i][j];
			v_new2[i][j]=yvel[i][j];
			v_new3[i][j]=yvel[i][j];
			v_new[i][j]=yvel[i][j];
			v_k1[i][j]=yvel[i][j];
			v_k2[i][j]=yvel[i][j];
			v_k3[i][j]=yvel[i][j];
			v_k[i][j]=yvel[i][j];
			v_n1[i][j]=yvel[i][j];
			v_n2[i][j]=yvel[i][j];
			v_n3[i][j]=yvel[i][j];
			v_n[i][j]=yvel[i][j];
			v_old1[i][j]=yvel[i][j];
			v_old2[i][j]=yvel[i][j];
			v_old3[i][j]=yvel[i][j];
			v_old[i][j]=yvel[i][j];
			p_new1[i][j]=Press[i][j];
			p_new2[i][j]=Press[i][j];
			p_new3[i][j]=Press[i][j];
			p_new[i][j]=Press[i][j];
			p_n[i][j]=Press[i][j];
		    dtau[i][j]=0.00005;
		    	
		}
	}
	double dj=1;

    for(double ti=0; ti<t1;ti+=dt)
    {
		
		//while( dj>=0.00001)
		for (int lop = 0; lop < 100; ++lop)
		{
		 	 

		 	/*Fourth order Runge-Kutta*/

		    RHS(N,JC,u_k,v_k,Press,zx,ey,ex,zy,Re,eps,ep,rcs,rus,rvs,rho1,rho2); 
		    
		    for (int i =1; i<N-1; i++)			// First step of RK
		    {    for (int j =1; j<N-1; j++)
		        {
		            RHSu[i][j]=(((-3*u_k[i][j]+4*u_n[i][j]-u_old[i][j])/(2*dt))+rus[i][j]);
		            
		    	    RHSv[i][j]=(((-3*v_k[i][j]+4*v_n[i][j]-v_old[i][j])/(2*dt))+rvs[i][j]);
		    	    
		            p_new1[i][j]=Press[i][j]+0.25*(dtau[i][j]*rcs[i][j]);
		            u_new1[i][j]=u_k[i][j]+0.25*(dtau[i][j]*RHSu[i][j]);
		            v_new1[i][j]=v_k[i][j]+0.25*(dtau[i][j]*RHSv[i][j]);
		        }
		    }
		    BC(N,u_new1,v_new1,p_new1);			// BC after first step RK
		    
		    RHS(N,JC,u_new1,v_new1,p_new1,zx,ey,ex,zy,Re,eps,ep,rcs,rus,rvs,rho1,rho2);  
		    
		    for (int i =1; i<N-1; i++)			// Second step of RK
		    {
		        for (int j =1; j<N-1; j++)
		        {
		            RHSu11[i][j]=(((-3*u_k1[i][j]+4*u_n1[i][j]-u_old1[i][j])/(2*dt))+rus[i][j]);
		            
		    	    RHSv11[i][j]=(((-3*v_k1[i][j]+4*v_n1[i][j]-v_old1[i][j])/(2*dt))+rvs[i][j]);
		    	    
		            p_new2[i][j]=Press[i][j]+0.33*(dtau[i][j]*rcs[i][j]);
		            u_new2[i][j]=u_k[i][j]+0.33*(dtau[i][j]*RHSu11[i][j]);
		            v_new2[i][j]=v_k[i][j]+0.33*(dtau[i][j]*RHSv11[i][j]);
		        }
		    }
		    BC(N,u_new2,v_new2,p_new2);			// BC after second step RK
		    
		    RHS(N,JC,u_new2,v_new2,p_new2,zx,ey,ex,zy,Re,eps,ep,rcs,rus,rvs,rho1,rho2); 
		    
		    for (int i =1; i<N-1; i++)			// Third step of RK
		    {
		        for (int j =1; j<N-1; j++)
		        {   
		            RHSu22[i][j]=(((-3*u_k2[i][j]+4*u_n2[i][j]-u_old2[i][j])/(2*dt))+rus[i][j]);
		            
		    	    RHSv22[i][j]=(((-3*v_k2[i][j]+4*v_n2[i][j]-v_old2[i][j])/(2*dt))+rvs[i][j]);
		    	    
		            p_new3[i][j]=Press[i][j]+0.5*(dtau[i][j]*rcs[i][j]);
		            u_new3[i][j]=u_k[i][j]+0.5*(dtau[i][j]*RHSu22[i][j]);
		            v_new3[i][j]=v_k[i][j]+0.5*(dtau[i][j]*RHSv22[i][j]);
		        }
		    }
		    BC(N,u_new3,v_new3,p_new3);			// BC after third step RK
		    
		    RHS(N,JC,u_new3,v_new3,p_new3,zx,ey,ex,zy,Re,eps,ep,rcs,rus,rvs,rho1,rho2); 
		    
		    for (int i =1; i<N-1; i++)			// Fourth step of RK
		    {
		        for (int j =1; j<N-1; j++)
		        {   
		            RHSu33[i][j]=(((-3*u_k3[i][j]+4*u_n3[i][j]-u_old3[i][j])/(2*dt))+rus[i][j]);
		            
		    	    RHSv33[i][j]=(((-3*v_k3[i][j]+4*v_n3[i][j]-v_old3[i][j])/(2*dt))+rvs[i][j]);
		            
		            p_new[i][j]=Press[i][j]+(dtau[i][j]*rcs[i][j]);
		            u_new[i][j]=u_k[i][j]+(dtau[i][j]*RHSu33[i][j]);
		            v_new[i][j]=v_k[i][j]+(dtau[i][j]*RHSv33[i][j]);
		        }
		    }
		    BC(N,u_new,v_new,p_new);			// BC after fourth step RK`

		   
		     for (int i =0; i<N; i++)			// Updating old values
		    {
		        for (int j =0; j<N; j++)
		        { 
				    u_old1[i][j]=u_n1[i][j];
				    u_n1[i][j]=u_k1[i][j];
				    u_k1[i][j]=u_new1[i][j];
				    u_old2[i][j]=u_n2[i][j];
				    u_n2[i][j]=u_k2[i][j];
				    u_k2[i][j]=u_new2[i][j];
				    u_old3[i][j]=u_n3[i][j];
				    u_n3[i][j]=u_k3[i][j];
				    u_k3[i][j]=u_new3[i][j];
				    u_old[i][j]=u_n[i][j];
				    u_n[i][j]=u_k[i][j];
				    u_k[i][j]=u_new[i][j];
				    v_old1[i][j]=v_n1[i][j];
				    v_n1[i][j]=v_k1[i][j];
				    v_k1[i][j]=v_new1[i][j];
				    v_old2[i][j]=v_n2[i][j];
				    v_n2[i][j]=v_k2[i][j];
				    v_k2[i][j]=v_new2[i][j];
				    v_old3[i][j]=v_n3[i][j];
				    v_n3[i][j]=v_k3[i][j];
				    v_k3[i][j]=v_new3[i][j];
				    v_old[i][j]=v_n[i][j];
				    v_n[i][j]=v_k[i][j];
				    v_k[i][j]=v_new[i][j];		    
				    Press[i][j]=p_new[i][j];
			  }
		    } 		   
		    
		    
		}  	
    }

   /* Writing Data to file */
	    ofstream xout("xvel.vtk");
		xout<<"# vtk DataFile Version 2.0\nVTK from matlab\n"<<"ASCII\n"<<"DATASET STRUCTURED_POINTS\n"<<"DIMENSIONS "<<N<<" "<<N<<" "<<1<<"\nSPACING 1 1 1 \n"<<"ORIGIN 0 0 0\n"<<"POINT_DATA "<<N*N<<"\nSCALARS xvel float 1\n"<<"LOOKUP_TABLE default\n";
		for (int i = 0; i < N; ++i)
		{
			for (int j = 0; j < N; ++j)
			{
				xout<<u_new[i][j]<<" ";
			}
		}
		
		xout.close();
		
		ofstream yout("yvel.vtk");
		yout<<"# vtk DataFile Version 2.0\nVTK from matlab\n"<<"ASCII\n"<<"DATASET STRUCTURED_POINTS\n"<<"DIMENSIONS "<<N<<" "<<N<<" "<<1<<"\nSPACING 1 1 1 \n"<<"ORIGIN 0 0 0\n"<<"POINT_DATA "<<N*N<<"\nSCALARS yvel float 1\n"<<"LOOKUP_TABLE default\n";
		for (int i = 0; i < N; ++i)
		{
			for (int j = 0; j < N; ++j)
			{
				yout<<yvel[i][j]<<" ";
			}
		}
		
		yout.close();
		ofstream pout("Press.vtk");
		pout<<"# vtk DataFile Version 2.0\nVTK from matlab\n"<<"ASCII\n"<<"DATASET STRUCTURED_POINTS\n"<<"DIMENSIONS "<<N<<" "<<N<<" "<<1<<"\nSPACING 1 1 1 \n"<<"ORIGIN 0 0 0\n"<<"POINT_DATA "<<N*N<<"\nSCALARS Pressure float 1\n"<<"LOOKUP_TABLE default\n";
		for (int i = 0; i < N; ++i)
		{
			for (int j = 0; j < N; ++j)
			{
				pout<<p_new[i][j]<<" ";
			}
		}
		
		pout.close();
	    return 0;
}
