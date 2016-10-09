/****Tyler Green Vortex****
***************************

***************************/
#include<iostream>
#include<fstream>
#include<math.h>
#include<stdlib.h>
#include"init.h"
#include"gridgen.h"
#include"metric.h"
#include"RHS.h"
#include"BC.h"

using namespace std;

int main()
{
	
	int N= 11;                          
	double st=1; 			  
	double Re=100;
	double eps=0.01;                         
	double L=30;
	double t=0;
	double t1=0.1;
	double Etd=0;
	double Es=0;
	double Es_st=0;
	double dt=0.005;
	double** x; double**y; double** xvel; double** yvel; double** xvel1; double** yvel1; double** Press;
	x=new double* [N];  y=new double* [N];  xvel=new double* [N];  yvel=new double* [N];  xvel1=new double* [N];  yvel1=new double* [N];  Press=new double* [N];		  	
	for (int i = 0; i < N; ++i)
	{
		x[i]=new double [N];  y[i]=new double [N];  xvel[i]=new double [N];  yvel[i]=new double [N];  xvel1[i]=new double [N];  yvel1[i]=new double [N];  Press[i]=new double [N];		  	
	}
	double norms, sum, norm_diff[N][N];    
	gridgen(N, st, x, y);  // Grid generation
	
	
	init( N, Re, t, t1, x, y, xvel, xvel1, yvel, yvel1, Press); 
	
	// Initial conditions
	
	
	double** u_new; double** u_new1; double** u_new2; double** u_new3; double** v_new; double** v_new1; double** v_new2; double** v_new3; double** p_new; double** p_new1; double** p_new2; double** p_new3;
	double** u_k; double** u_k1; double** u_k2; double** u_k3; double** v_k; double** v_k1; double** v_k2; double** v_k3; double** p_k; double** p_k1; double** p_k2; double** p_k3;
	double** u_old; double** u_old1; double** u_old2; double** u_old3; double** v_old; double** v_old1; double** v_old2; double** v_old3; double** u_n; double** u_n1; double** u_n2; double** u_n3;
	double** v_n; double** v_n1; double** v_n2; double** v_n3; double** p_n; double** dtau; double** JC; double** ex; double** ey; double** zx; double** zy; double** rus; double** rvs; double** rcs;
	double** RHSu11; double** RHSu22; double** RHSu33; double** RHSv11; double** RHSv22; double** RHSv33; double** p1; double** p2; double** p3; double** RHSu; double** RHSv; double** rho1; double** rho2;

	u_new=new double* [N];  u_new1=new double* [N];  u_new2=new double* [N];  u_new3=new double* [N];  v_new=new double* [N];  v_new1=new double* [N];  v_new2=new double* [N];  v_new3=new double* [N];  p_new=new double* [N];  p_new1=new double* [N];  p_new2=new double* [N];  p_new3=new double* [N];
	u_k=new double* [N];  u_k1=new double* [N];  u_k2=new double* [N];  u_k3=new double* [N];  v_k=new double* [N];  v_k1=new double* [N];  v_k2=new double* [N];  v_k3=new double* [N];  p_k=new double* [N];  p_k1=new double* [N];  p_k2=new double* [N];  p_k3=new double* [N];
	u_old=new double* [N];  u_old1=new double* [N];  u_old2=new double* [N];  u_old3=new double* [N];  v_old=new double* [N];  v_old1=new double* [N];  v_old2=new double* [N];  v_old3=new double* [N];  u_n=new double* [N];  u_n1=new double* [N];  u_n2=new double* [N];  u_n3=new double* [N];
	v_n=new double* [N];  v_n1=new double* [N];  v_n2=new double* [N];  v_n3=new double* [N];  p_n=new double* [N];  dtau=new double* [N];  JC=new double* [N];  ex=new double* [N];  ey=new double* [N];  zx=new double* [N];  zy=new double* [N];  rus=new double* [N];  rvs=new double* [N];  rcs=new double* [N];
	RHSu11=new double* [N];  RHSu22=new double* [N];  RHSu33=new double* [N];  RHSv11=new double* [N];  RHSv22=new double* [N];  RHSv33=new double* [N];  p1=new double* [N];  p2=new double* [N];  p3=new double* [N];  RHSu=new double* [N];  RHSv=new double* [N];rho1=new double* [N];rho2=new double* [N];

	for (int i = 0; i < N; ++i)
	{
		 u_new[i]=new double [N];  u_new1[i]=new double [N];  u_new2[i]=new double [N];  u_new3[i]=new double [N];  v_new[i]=new double [N];  v_new1[i]=new double [N];  v_new2[i]=new double [N];  v_new3[i]=new double [N];  p_new[i]=new double [N];  p_new1[i]=new double [N];  p_new2[i]=new double [N];  p_new3[i]=new double [N];
		 u_k[i]=new double [N];  u_k1[i]=new double [N];  u_k2[i]=new double [N];  u_k3[i]=new double [N];  v_k[i]=new double [N];  v_k1[i]=new double [N];  v_k2[i]=new double [N];  v_k3[i]=new double [N];  p_k[i]=new double [N];  p_k1[i]=new double [N];  p_k2[i]=new double [N];  p_k3[i]=new double [N];
		 u_old[i]=new double [N];  u_old1[i]=new double [N];  u_old2[i]=new double [N];  u_old3[i]=new double [N];  v_old[i]=new double [N];  v_old1[i]=new double [N];  v_old2[i]=new double [N];  v_old3[i]=new double [N];  u_n[i]=new double [N];  u_n1[i]=new double [N];  u_n2[i]=new double [N];  u_n3[i]=new double [N];
		 v_n[i]=new double [N];  v_n1[i]=new double [N];  v_n2[i]=new double [N];  v_n3[i]=new double [N];  p_n[i]=new double [N];  dtau[i]=new double [N];  JC[i]=new double [N];  ex[i]=new double [N];  ey[i]=new double [N];  zx[i]=new double [N];  zy[i]=new double [N];  rus[i]=new double [N];  rvs[i]=new double [N];  rcs[i]=new double [N];
		 RHSu11[i]=new double [N];  RHSu22[i]=new double [N];  RHSu33[i]=new double [N];  RHSv11[i]=new double [N];  RHSv22[i]=new double [N];  RHSv33[i]=new double [N];  p1[i]=new double [N];  p2[i]=new double [N];  p3[i]=new double [N];  RHSu[i]=new double [N];  RHSv[i]=new double [N];rho1[i]=new double [N];rho2[i]=new double [N];
	}

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
		 	sum=0;   

		 	/*Fourth order Runge-Kutta*/

		    RHS(N,JC,u_k,v_k,Press,zx,ey,ex,zy,Re,eps,rcs,rus,rvs,rho1,rho2); 
		    
		    for (int i =1; i<N-1; i++)			// First step of RK
		    {    for (int j =1; j<N-1; j++)
		        {
		            RHSu[i][j]=(((-3*u_k[i][j]+4*u_n[i][j]-u_old[i][j])/(2*dt))+rus[i][j]);
		            
		    	    RHSv[i][j]=(((-3*v_k[i][j]+4*v_n[i][j]-v_old[i][j])/(2*dt))+rvs[i][j]);
		    	    
		            p_new[i][j]=Press[i][j]+0.25*(dtau[i][j]*rcs[i][j]);
		            u_new[i][j]=u_k[i][j]+0.25*(dtau[i][j]*RHSu[i][j]);
		            v_new[i][j]=v_k[i][j]+0.25*(dtau[i][j]*RHSv[i][j]);
		        }
		    }
		    BC(N,u_new,v_new,p_new);			// BC after first step RK
		    
		    RHS(N,JC,u_new,v_new,p_new,zx,ey,ex,zy,Re,eps,rcs,rus,rvs,rho1,rho2); 
		    
		    for (int i =1; i<N-1; i++)			// Second step of RK
		    {
		        for (int j =1; j<N-1; j++)
		        {
		            RHSu11[i][j]=(((-3*u_k1[i][j]+4*u_n1[i][j]-u_old1[i][j])/(2*dt))+rus[i][j]);
		            
		    	    RHSv11[i][j]=(((-3*v_k1[i][j]+4*v_n1[i][j]-v_old1[i][j])/(2*dt))+rvs[i][j]);
		    	    
		            p_new[i][j]=Press[i][j]+0.33*(dtau[i][j]*rcs[i][j]);
		            u_new[i][j]=u_k[i][j]+0.33*(dtau[i][j]*RHSu11[i][j]);
		            v_new[i][j]=v_k[i][j]+0.33*(dtau[i][j]*RHSv11[i][j]);
		        }
		    }
		    BC(N,u_new,v_new,p_new);			// BC after second step RK
		    
		    RHS(N,JC,u_new,v_new,p_new,zx,ey,ex,zy,Re,eps,rcs,rus,rvs,rho1,rho2); 
		    
		    for (int i =1; i<N-1; i++)			// Third step of RK
		    {
		        for (int j =1; j<N-1; j++)
		        {   
		            RHSu22[i][j]=(((-3*u_k2[i][j]+4*u_n2[i][j]-u_old2[i][j])/(2*dt))+rus[i][j]);
		            
		    	    RHSv22[i][j]=(((-3*v_k2[i][j]+4*v_n2[i][j]-v_old2[i][j])/(2*dt))+rvs[i][j]);
		    	    
		            p_new[i][j]=Press[i][j]+0.5*(dtau[i][j]*rcs[i][j]);
		            u_new[i][j]=u_k[i][j]+0.5*(dtau[i][j]*RHSu22[i][j]);
		            v_new[i][j]=v_k[i][j]+0.5*(dtau[i][j]*RHSv22[i][j]);
		        }
		    }
		    BC(N,u_new,v_new,p_new);			// BC after third step RK
		    
		    RHS(N,JC,u_new3,v_new3,p_new3,zx,ey,ex,zy,Re,eps,rcs,rus,rvs,rho1,rho2); 
		    
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
				    p1[i][j]=p_new1[i][j];
				    p2[i][j]=p_new2[i][j];
				    p3[i][j]=p_new3[i][j];
				    Press[i][j]=p_new[i][j];
			  }
		    } 		   
		    
		    
		}  	
    }

    /* Writing Data to file */
	    ofstream xout("xvel.csv");
		for (int i = 0; i < N; ++i)
		{
			for (int j = 0; j < N; ++j)
			{
				xout<<u_new[i][j]<<",";
			}
			xout<<endl;
		}
		
		xout.close();
		ofstream yout("yvel.csv");
		for (int i = 0; i < N; ++i)
		{
			for (int j = 0; j < N; ++j)
			{
				yout<<v_new[i][j]<<",";
			}
			yout<<endl;
		}
		
		yout.close();
		ofstream pout("Press.csv");
		for (int i = 0; i < N; ++i)
		{
			for (int j = 0; j < N; ++j)
			{
				pout<<p_new[i][j]<<",";
			}
			pout<<endl;
		}
		
		pout.close();
	    return 0;
}
