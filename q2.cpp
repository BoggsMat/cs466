#include <iostream>
#include <math.h> 

using namespace std; 

#define n 4

void printCol(double vec[n])
	{
	for (int i = 0; i < n; ++i)
		{
		cout << vec[i] << ' '; 	
		}
	cout << endl; 
	}

void printMatrix(double A[n][n])
	{
    for(int i=0;i<n;i++)
    	{
        for(int j=0;j<n;j++)
            cout << A[i][j] << ' '; 
        cout << endl;
    	}
	}

void copyColumn(double newCol[n], double orig[n][n], int colNum)
	{
	for(int i = 0; i < n; ++i)
		newCol[i] =orig[i][colNum]; 
	}

void copyFromt(double newCol[n], double orig[n])
	{
	for(int i = 0; i < n; ++i)
		newCol[i] =orig[i]; 
	}

double dotprod(double x[n],double y[n])
	{
    double r=0.0;
    for(int k=0;k<n;k++) 
    	r+=x[k]*y[k];

    return r;
	}

double calcNorm(double vec[n])
	{
	double prod = dotprod(vec, vec);

	return sqrt(prod);
	}

void subtractVectors(double vec[n], double vec2[n], double result[n])
	{
	for(int i = 0; i < n; ++i)
		result[i] = vec[i] - vec2[i]; 
	}

void multiVector(double vec[n], double temp[n], double val)
	{
	for (int i = 0; i < n; ++i)
		temp[i] = vec[i] * val; 
	}

void divideVector(double vec[n], double val)
	{
	for (int i = 0; i < n; ++i)
		{
		vec[i] /= val; 
		}
	}

void createQ(double Q[n][n], double v[n][n])
	{
	for(int i = 0; i < n; ++i)
		for(int j = 0; j < n; ++j)
			{
			Q[i][j] = v[j][i]; 
			}
	}

void gramSchmidt(double A[n][n], double Q[n][n], double R[n][n])
	{
	double t[n][n][n], v[n][n], temp[n]; 
	int i = 0, j = 0; 
	double norm;  

	for(i = 0; i < n; ++i)
		{	
		for(j = 0; j <= i; ++j)
			{
			// Find ti,j
			if (j == 0)
				copyColumn(t[0][i], A, i);
			else
				{
				multiVector(v[j-1], temp, dotprod(v[j-1], t[j-1][i]));
				subtractVectors(t[j-1][i], temp, t[j][i]);

				}
			}
		// Calculate vn
		copyFromt(v[i], t[i][i]);
		norm = calcNorm(v[i]); 
		divideVector(v[i], norm); 
		}

	createQ(Q, v); 

	for (i = 0; i < n; ++i)
		for(j = 0; j < n; ++j)
			{
			if (i < j)
				R[i][j] = dotprod(v[i], t[i][j]); 
			else if(i == j)
				R[i][j] = calcNorm(t[i][i]); 
			else 
				R[i][j] = 0; 
			}

	}

void hilbert(double H[n][n])
	{
	for (int i = 1; i <= n; i++)
		for (int j = 1; j <= n; j++)
			H[i-1][j-1] = 1.0/(i+j-1.0); 
	}

double fnorm(double A[n][n])
	{
	double r=0.0;
	for(int i=0;i<n;i++)
		for(int j=0;j<n;j++)
			r+=A[i][j]*A[i][j];
		
	
	return sqrt(r);
	}

int main()
	{
	double Q[n][n], R[n][n], H[n][n], T[n][n]; 

	// Create Hilbert Matrix 
	hilbert(H); 

	// Calculate the QR factorization 
	gramSchmidt(H, Q, R); 

	cout << "H:" << endl; 
	printMatrix(H); 

	cout << endl << "Q:" << endl; 
	printMatrix(Q); 

	cout << endl << "R:" << endl; 
	printMatrix(R); 

	return 0;
	}