#include <iostream> 
#include <math.h> 

using namespace std; 

double y(double t)
	{
	return t*t*(exp(t) - exp(1));
	}

double f(double t, double y)
	{
	return (2/t)*y + t*t*exp(t); 
	}

double df(double t, double y)
	{
	return -2*(y/(t*t)) + (2/t)*f(t,y) + 2*t*exp(t) + t*t*exp(t); 
	}

double d2f(double t, double y)
	{
	return (4/t*t*t)*y - (4/t*t)*f(t,y) + exp(t)*t*t + (2/t)*df(t,y) + 4*exp(t)*t + 2*exp(t);
	}

double d3f(double t, double y)
	{
	return (-12/t*t*t*t)*y +(12/t*t*t)*f(t,y) - (6/t*t)*df(t,y) + (2/t)*d2f(t,y) + exp(t)*t*t + 6*exp(t)*t + 6*exp(t);
	}

double linearInterp(double xin, double x1, double y1, double x2, double y2)
	{
	double slope = (y2 - y1)/(x2 -x1); 

	return slope*(xin - x1) + y1; 
	}

double hermite(double xin)
	{
	int n = 3; 
	double z[2*n+1];
	double Q[2*n+1][2*n+1]; 
	double x[] = {1, 1.5, 2}; 
	double fun[] = {f(1, 1), f(1.5, 1), f(2, 1)}; 
	double dfun[] = {df(1,1), df(1.5, 1) , df(2, 1)}; 

	// Step 1
	for (int i = 0; i < n; ++i)
		{
		// Step 2
		z[2*i] = z[2*i+1] = x[i]; 
		Q[2*i][0] = Q[2*i+1][0]= fun[i]; 
		Q[2*i+1][1] = dfun[i]; 

		// Step 3
		if (i != 0)
			Q[2*i][1] = (Q[2*i][0] - Q[2*i-1][0])/(z[2*i] - z[2*i-1]); 
		}

	// Step 4
	for (int j = 2; j < (2*n+1); ++j)
		for (int k = 2; k < n; ++k)
			{
			Q[j][k] = (Q[j][k-1] - Q[j-1][k-1])/(z[j] - z[j-k]); 
			}

	return Q[0][0] + Q[1][1]*xin + Q[2][2]*(xin*xin); 
	}

void eulersMethod(double a, double h, double alpha)
	{
	double t = a; 
	double w = alpha; 
	double actual = y(t);

	for(int i = 1; i <= 10; ++i)
		{
		w += h*f(t, w); 
		t += h; 
		actual = y(t);
		cout << t << "  "  << w << "  "  << actual << endl; 
		}
	}

void taylorsMethod2(double a, double h, double alpha)
	{
	double t = a; 
	double w = alpha; 
	double actual = y(t); 

	for (int i = 1; i <= 10; ++i)
		{
		w += h*(f(t, w) + (h/2)*df(t,w));
		t += h;  
		double actual = y(t); 
		cout << t << "  "  << w << "  "  << actual << endl; 
		}
	}

void taylorsMethod4(double a, double h, double alpha)
	{
	double t = a; 
	double w = alpha; 
	double actual = y(t); 

	for (int i = 1; i <= 10; ++i)
		{
		w += h*(f(t, w) + (h/2)*df(t,w) + (h*h/6)*d2f(t,w) + (h*h*h/24)*d3f(t,w));
		t += h;  
		double actual = y(t); 
		cout << t << "  "  << w << "  "  << actual << endl; 
		}
	}

int main()
	{
	// 5.2
		// 9a
		cout << "Eulers Method:" << endl; 
		eulersMethod(1, 0.1, 0);
		cout << endl; 

		// 9b
		cout << "Linear Interpolation" << endl; 
		cout << "f(1.04) = " << linearInterp(1.04, 1, 0, 1.1, 0.2718) << endl;
		cout << "f(1.55) = " << linearInterp(1.55, 1.5, 3.1875, 1.6, 4.6208) << endl; 
		cout << "f(1.97) = " << linearInterp(1.97, 1.9, 11.748, 2, 15.3982) << endl << endl; 

	//5.3 
		// 9a 
		cout << "Taylors Method order two" << endl; 
		taylorsMethod2(1, 0.1, 0);
		cout << endl; 

		// 9b
		cout << "Linear Interpolation" << endl; 
		cout << "la(1.04) = " << linearInterp(1.04, 1, 0, 1.1, 0.3398) << endl;
		cout << "la(1.55) = " << linearInterp(1.55, 1.5, 3.911, 1.6, 5.6431) << endl; 
		cout << "la(1.97) = " << linearInterp(1.97, 1.9, 14.1527, 2, 18.47) << endl << endl; 

		// 9c 
		cout << "Taylors Method order four" << endl; 
		taylorsMethod4(1, 0.1, 0);
		cout << endl; 

		// 9d 
		cout << "Linear Interpolation" << endl; 
		cout << "h(1.04) = " << hermite(1.04) << endl;
		cout << "f(1.55) = " << linearInterp(1.55, 1.5, 3.9618, 1.6, 5.7116) << endl; 
		cout << "f(1.97) = " << linearInterp(1.97, 1.9, 14.299, 2, 18.6536) << endl << endl; 
	return 0; 
	}