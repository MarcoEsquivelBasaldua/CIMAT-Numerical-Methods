/*
*   Code name: Power_Iteration.c
*   Author: Marco Antonio Esquivel Basaldua
*   Date: 05/09/19
*
*   Code description:
*   This Code captures Matrix A and finds its maximum eigenvalue and its eigenvector assotiated using the power iteration method
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>


void matXvec(double**,double*,double*,int);                                 // Function to multiply a matrix times a vector
void normV(double*,int);                                                   // This function calculates the norm of the given vector
double Power_Iteration(double**,double*,double*,int,double,int);          // Power iteration method to find maximum eigenvalue and eigenvector assotiated


int main(int argc , char* argv[] ){
    int n,m;        // Size of matrix A
    double **A;     // Matrix A
    int i,j;        // Integers used in all for cycles

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

    double* V0 = malloc(sizeof(double) * n);        // VO is nedded to find the corresponding eigenvector
    double* V1 = malloc(sizeof(double) * n);        // Corresponding eigenvector will be stored here

    double lamda = 0;                        // Maximum eigenvalue will be stored here
    double toler = 1e-5;                    // A tolerance will determine if the code has found the solutions
    int maxIter = 10000;                   // This value is the maximum of iterations the code will run. This avoids the code to execute forever

    lamda = Power_Iteration(A,V0,V1,n,toler,maxIter);


    // Memory is freed
    for(int i=0;i<n;i++) free(A[i]);
    free(A);
    free(V0);
    free(V1);

    return 0;
}

// This function calculates the norm of the given vector
void normV(double *V,int n){
/*
*   Inputs:
*       - double *V: vector to be normalized, its corresponding vector normalized will be passed in this very same location
*       - int n: number of rows in vector V
*
*   Outputs:
*       - None
*/
    double den = 0;
    for(int i=0; i<n; i++) den += V[i]*V[i];

    den = sqrt(den);

    for(int i=0; i<n; i++) V[i] = V[i]/den;
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

double Power_Iteration(double **A, double *V0, double *V1, int n, double toler, int maxIter){
/*
*   Intputs:
*       - double **A: pointer to matrix whom maximum eigenvalue will be here calculated
*       - double *VO: pointer to vector
*       - double *V1: pointer to vector where eigenvector will be passed by deferring
*       - int n: number of rows in matrix A and vector VO and V1
*       - double toler: tolerance near to cero to determine if the code has found the solutions
*       - int max_iter: This value is the maximum of iterations the code will run. This avoids the code to execute forever
*
*   Outputs:
*       - double lamda: maximum eigenvalue of matrix A
*/
    int i, j;
    double lambda = 0;
    double lambda_old = 0;

    srand(time(NULL));
    // Initialize V0 with random numbers
    for(i=0; i<n ; i++) V0[i] = rand()%100+1;
    //V0[0] = 1;


    // Iterations start to find maximum eigenvalue and eigenvector
    for(int iter=1; iter<maxIter; iter++){
        normV(V0,n);            // Normalize V0 to avoid working with large values
        matXvec(A,V0,V1,n);     // We get V1 by multiplying A*V0

        double num = 0;
        double den = 0;

        // The value proposed for the eigenvalue is proposed
        for(i=0; i<n; i++){
            num += V1[i]*V1[i];
            den += V1[i]*V0[i];
        }

        lambda = num/den;

        // lambda will be set as the eigenvalue only if the difference with the previous value is far too small,
        // this would mean that the value has converged
        if(fabs(lambda-lambda_old)<toler){
            printf("Iterations taken: %d\n", iter);     // This will help know how many iterations it took to get to the results
            normV(V1,n);                               // Before leaving this loop the eigenvector found is normalized
            break;                                    // break to leave this foor loop
        }
        else{                                   // In the case the difference between the current proposal to eigenvalue and the last one is significant
            lambda_old = lambda;               // lambda_old will take the value of lambda to start next iteration
            for(i=0;i<n;i++) V0[i] = V1[i];   // as VO will take values of V1
        }
    }

    // Results are displayed
    printf("Maximum eigenvalue is:\n %lf\n", lambda);

    printf("Eigenvector assotiated to maximum eigenvalue is:\n");
    for(i=0; i<n; i++) printf("%.10lf\n",V1[i]);

    return lambda;
}
