#include<math.h>

using namespace std;
// Initializing U, V,  Pressure
void init(int N,  double Re, double t, double t1, vector<vector<double> > &x, vector<vector<double> > &y, vector<vector<double> > &xvel, vector<vector<double> > &xvel1, vector<vector<double> > &yvel, vector<vector<double> > &yvel1, vector<vector<double> > &Press)
{
	double Press1[N][N];
	for (int i=0; i<N; i++)
	{
	    for (int j=0; j<N; j++)
    	    {
	        xvel[i][j]=-exp((-2*t/Re))*cos(x[1][i])*sin(y[j][1]);
	        yvel[i][j]=exp((-2*t)/Re)*cos(y[j][1])*sin(x[1][i]);
	        Press[i][j]=-(((cos(2*x[1][i]))+cos(2*y[j][1]))*exp((-4*t)/Re))/4;
	        xvel1[i][j]=-exp((-2*t1)/Re)*cos(x[1][i])*sin(y[j][1]);
	        yvel1[i][j]=exp((-2*t1)/Re)*cos(y[j][1])*sin(x[1][i]);
	        Press1[i][j]=-((cos(2*x[1][i])+cos(2*y[j][1]))*exp((-4*t1)/Re))/4;
		        
	    }
	}
	
}
