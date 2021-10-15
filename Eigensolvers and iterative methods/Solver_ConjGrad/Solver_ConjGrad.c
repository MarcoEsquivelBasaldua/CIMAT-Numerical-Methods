/*
*   Code name: Solver_ConjGrad.c
*   Author: Marco Antonio Esquivel Basaldua
*   Date: 21/09/19
*
*   Code description:
*   This Code captures Matrix A,
*   captures the values of solutions to the n linear equations in the matrix b,
*   If A is a square Matrix and b has sufficient values as solutions for the linear system:
*      calculates the values of x from the equation Ax = b, by the conjugate gradient method that is an iterative method
*      used to find the values of the unknowns
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

void Solver_ConjGrad(double**,double**,int,double,int);  // Function to find the values of the unknowns by Jacobi's method
double dotProd(double *V1,double *V2,int n);            // Calculates the dot product given two vectors
void matXvec(double**,double*,double*,int);            // Function to multiply a matrix times a vector
double normV(double*,int);                            // Finds the euclidean norm of a vector

int main( int argc , char* argv[] ){
    int n,m;                // Size of matrix A
    int k,l;                // Size of matrix b (solutions to linear equations)
    double **A;             // Matrix A
    double **b;             // Values of solutions to the n linear equations
    double toler = 1e-6;    // A tolerance will determine if the code has found the solutions
    int max_iter = 10000;   // This value is the maximum of iterations the code will run. This avoids the code to execute forever
    int i,j;                // Integers used in all for cycles

    FILE* fin = NULL;
    // File with matrix A is read and displayed
	fin = fopen( argv[ 1 ] , "r" );
	if(  !fin  ){
		printf("Error: No se abrio %s\n" , argv[ 1 ] );
	}

    fscanf(fin, "%d %d", &n, &m );
    printf("Matrix A:\n");
	printf( "Rows: %d, Columns: %d\n" , n, m );
    A = malloc(n * sizeof *A);
    for (i=0; i<n; i++) A[i] = malloc(m * sizeof *A[i]);

    for(i=0; i<n; i++) for(j=0; j<m; j++) fscanf(fin, "%lf", &A[i][j]);

    if( n<=20 && m<=20){                                             // Matrix A is only displayed if it is smaller than 20 rows and cols
        for(i=0; i<n; i++){
		    for(j=0; j<m; j++){
			    printf("%lf ", A[i][j]);
		    }
		    printf("\n");
	    }
    }
    else printf("Matrix A is too large to be displayed\n");
    printf("\n");
	fclose( fin );

    // File with matrix b is read and displayed
	fin = fopen( argv[ 2 ] , "r" );
	if(  !fin  ){
		printf( "Error: No se abrio %s\n" , argv[ 1 ] );
	}

    fscanf( fin , "%d %d", &k , &l );
    printf("Matrix b:\n");
	printf( "Rows: %d, Columns: %d\n" , k, l );
    b = malloc(k * sizeof *b);
    for (i=0; i<k; i++) b[i] = malloc(l * sizeof *b[i]);

    for(i=0; i<k; i++) for(j=0; j<l; j++) fscanf(fin, "%lf", &b[i][j]);

    if( k<=20 && l<=20){                                             // Matrix b is only displayed if it is smaller than 20 rows and cols
        for(i=0; i<k; i++){
		    for(j=0; j<l; j++){
			    printf("%lf ", b[i][j]);
		    }
		    printf("\n");
	    }
    }
    else printf("Matrix b is too large to be displayed\n");
    printf("\n");
	fclose( fin );

    // Time excecution starts
    clock_t t_ini, t_fin;
    double secs;
    t_ini = clock();

    Solver_ConjGrad(A,b,n,toler,max_iter);
    printf("\n");

    // Time excecution is displayed
    t_fin = clock();
    secs = (double)(t_fin - t_ini) / CLOCKS_PER_SEC;
    printf("%.16g miliseconds taken in the algorithm\n", secs * 1000.0);

    // Memory is freed
    for(int i=0;i<n;i++) free(A[i]);
    free(A);
    for(int i=0;i<n;i++) free(b[i]);
    free(b);
    return 0;
}

void Solver_ConjGrad(double** A, double** b,int n, double toler, int maxIter){
/*
*   Inputs:
*       - double **A: Pointer to the matrix A
*       - double **b: Pointer to the matrix b
*       - int n: number of rows in A
*       - double toler: tolerance near to cero to determine if the code has found the solutions
*       - int max_iter: This value is the maximum of iterations the code will run. This avoids the code to execute forever
*/
    double* x = calloc(sizeof(double) , n);
    double* rold = malloc(sizeof(double) * n);
    double* omega = malloc(sizeof(double) * n);
    double* pold = malloc(sizeof(double) * n);

    for(int i=0; i<n; i++){
        rold[i] = b[i][0];
        pold[i] = b[i][0];
    }

    double alpha, betha;

    for(int iter=1; iter<=maxIter; iter++){

        matXvec(A,pold,omega,n);
        double num,den;
        num = dotProd(pold,rold,n);
        den = dotProd(pold,omega,n);

        alpha = num/den;

        for(int i=0; i<n; i++){
            x[i] = x[i] + alpha*pold[i];
            rold[i] = rold[i] - alpha*omega[i];
        }

        double norm = normV(rold,n);
        //printf("%lf\n",norm);
        if(norm <= toler){
            printf("Iterations taken: %d\n",iter);
            break;
        }
        else{
            num = dotProd(pold,rold,n);
            den = dotProd(pold,pold,n);
            betha = num/den;

            for(int i=0; i<n; i++) pold[i] = rold[i] + betha*pold[i];
        }
    }

    // Solutions to unknowns are displayed
    printf("The values of the unknowns are:\n");
    for(int i=0; i<n; i++) printf("%.10lf\n",x[i]);

    free(x);
    free(rold);
    free(pold);
    free(omega);
}

double dotProd(double *V1,double *V2,int n){
/*
*   Inputs:
*       - double *V1: first vector in the dot product operation
*       - double *V2: second vector in the dot product operation
*       - int n: dimension of vectors V1 and V2
*
*   Output:
*       - double: result on the dot product operation V1*V2
*/
    double Prod=0;

    for(int i=0;i<n;i++) Prod += V1[i]*V2[i];

    return Prod;
}

// Function to multiply a matrix times a vector
void matXvec(double **A, double *V0, double* V1, int n){
/*
*   Inputs:
*       - double **A: pointer to matrix to be multiplied
*       - double *VO: pointer to vector to be multiplied
*       - double *V1: multiplication result will be passed here by deferring
*       - int n: number of rows in matrix A and vector VO and V1
*
*   Outputs:
*       - None
*/
    for(int i=0; i<n; i++){
        V1[i] = 0;
        for(int j=0;j<n;j++){
            V1[i] += A[i][j]*V0[j];
        }
    }
}

// Finds the euclidean norm of a vector
double normV(double* V,int n){
/*
*   Intputs:
*       - double* V: vector to calculate it norm
*       - int n: number of entries in vector V
*
*   Outputs:
*       - double : euclidean norm of vector V
*/
    double norm = 0;
    for(int i=0; i<n; i++) norm += V[i]*V[i];

    norm = sqrt(norm);

    return norm;
}
