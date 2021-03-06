#include <iostream>
#include <math.h> 

using namespace std; 

#define n 6

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

void multiVector(double vec[n], double result[n], double val)
	{
	for (int i = 0; i < n; ++i)
		result[i] = vec[i] * val; 
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

double fnorm(double A[n][n])
	{
	double r=0.0;
	for(int i=0;i<n;i++)
		for(int j=0;j<n;j++)
			r+=A[i][j]*A[i][j];
		
	
	return sqrt(r);
	}

void multiMatrix(double A[n][n], double B[n][n], double result[n][n])
	{
	double temp = 0.0; 

	for (int k = 0; k < n; ++k)
		{
		for (int i = 0; i < n; ++i)
			{
			for (int j = 0; j < n; ++j)
				{
				temp += A[k][j]*B[j][i];
				}
			result[i][k] = temp; 
			temp = 0.0; 
			}
		}
	}

void subtractMatrix(double A[n][n], double B[n][n], double result[n][n])
	{
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j)
			result[i][j] = A[i][j] - B[i][j]; 
	}

void AtA(double B[n][n], double A[n][n])
	{
	bzero(B,sizeof(double)*n*n);
    for(int k=0;k<n;k++)
        for(int i=0;i<n;i++)
            for(int j=0;j<n;j++) 
                B[i][j] += A[k][i]*A[k][j];
	}

void AAt(double B[n][n], double A[n][n])
	{
	bzero(B,sizeof(double)*n*n);
    for(int k=0;k<n;k++)
        for(int i=0;i<n;i++)
            for(int j=0;j<n;j++) 
                B[k][i] += A[k][j]*A[i][j];
	}

int main()
	{
	double Q[n][n], R[n][n], H[n][n], T[n][n], A[n][n], I[n][n]; 

	// Create Identity 
	for (int i = 0; i < n; ++i)
		for(int j = 0; j < n; ++j)
			i == j ? I[i][j] = 1 : I[i][j] = 0; 

	// Create Hilbert Matrix 
	for (int i = 1; i <= n; i++)
		for (int j = 1; j <= n; j++)
			H[i-1][j-1] = 1.0/(i+j-1.0); 

	// Calculate the QR factorization 
	gramSchmidt(H, Q, R); 

	cout << "H:" << endl; 
	printMatrix(H); 

	cout << endl << "Q:" << endl; 
	printMatrix(Q); 

	cout << endl << "R:" << endl; 
	printMatrix(R); 

	// ||H-QR||
	cout << endl << "||H - QR|| " << endl; 
	multiMatrix(Q,R,A);
	subtractMatrix(H,A,T);
	cout << fnorm(T) << endl; 

	// ||QtQ - I||
	cout << endl << "||QtQ - I||" << endl; 
	AtA(T, Q); 
	subtractMatrix(T, I, A); 
	cout << fnorm(A) << endl; 

	// ||QQt - I||
	cout << endl << "||QQt - I||" << endl; 
	AAt(T, Q); 
	subtractMatrix(T, I, A); 
	cout << fnorm(A) << endl; 


	return 0;
	}