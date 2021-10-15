/*
*   Code name: Power_Iteration_kBiggest.c
*   Author: Marco Antonio Esquivel Basaldua
*   Date: 14/09/19
*
*   Code description:
*   This Code captures Matrix A and finds its k maximum eigenvalues and their
*   eigenvectors assotiated using the power iteration method.
*
*   The results are stored in two .txt files:
*       "Eigenvalues.txt"
*       "Eigenvectors.txt"
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define K 15    // This is the number of eigenvalues we are interested in

void matXvec(double**,double*,double*,int);                              // Function to multiply a matrix times a vector
void Power_Iteration(double**,double*,double**,int,double,int,int);     // Power Iteration method to find k biggest eigenvalues and eigenvectors assotiated
void normV(double*,int);                                               // Finds the euclidean norm from a vector
double dotProd(double*,double*,int);                                  // Calculates the dot product given two vectors

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

    double toler = 1e-12;                    // A tolerance will determine if the code has found the solutions
    int maxIter = 10000;                   // This value is the maximum of iterations the code will run. This avoids the code to execute forever

    double* lamda = malloc(sizeof(double) * K);        // All eigenvalues calculated will be here stored
    double **FI;                                       // Eigenvectors calculated will be here stored
    FI = malloc(K * sizeof *FI);
    for (int i=0; i<K; i++) FI[i] = malloc(n * sizeof *FI[i]);

    // Time excecution starts
    clock_t t_ini, t_fin;
    double secs;
    t_ini = clock();

    Power_Iteration(A,lamda,FI,n,toler,maxIter,K);
    // Time excecution is displayed
    t_fin = clock();
    secs = (double)(t_fin - t_ini) / CLOCKS_PER_SEC;
    printf("%.16g miliseconds taken in the algorithm\n", secs * 1000.0);

    // Results are displayed
    printf("File with Eigenvalues is generated.\n\n");
    /////////
    char *nombre = "Eigenvalues.txt";     // Nombre que se le dará a la imagen generada

    FILE* file = NULL;
	file = fopen( nombre, "w" );
	if(  !file  ){
		printf("Error: No se abrio %s\n" , nombre );
	}

    fprintf(file,"Minimun eigenvalues are:\n");
    for(int i=0; i<K; i++) fprintf(file,"%lf\n", lamda[i]);
    fclose(file);
    ////

    printf("File with Eigenvectors is generated.\n");
    /////////
    char *nombre2 = "Eigenvectors.txt";     // Nombre que se le dará a la imagen generada

    FILE* file2 = NULL;
	file2 = fopen( nombre2, "w" );
	if(  !file2  ){
		printf("Error: No se abrio %s\n" , nombre2 );
	}

    fprintf(file2,"Eigenvectors assotiated:\n");
    for(int j=0;j<n;j++){
        for(int i=0;i<K;i++) fprintf(file2,"%lf  ",FI[i][j]);
        fprintf(file2,"\n");
    }
    fclose(file2);
    ////


    // Memory is freed
    for(int i=0;i<n;i++) free(A[i]);
    free(A);
    for(int i=0;i<K;i++) free(FI[i]);
    free(FI);
    free(lamda);

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

    //for(int i =0;i<n;i++) printf("%lf\n",V[i]);
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

void Power_Iteration(double** A,double* lamda,double** FI,int n,double toler,int maxIter,int k){
/*
*   Inputs:
*       - double **A: square matrix to calculate it lowest eigenvalue and assotiated eigen vector
*       - double *lamda: calculated k eigenvalues will be stored here as a vector
*       - double **FI: mxn matrix where eigenvectors will be stored by rows
*       - int n: dimension of square matrix A
*       - double toler: tolerance near to cero to determine if the code has found the solutions
*       - int max_iter: This value is the maximum of iterations the code will run. This avoids the code to execute forever
*       - int k: given amount of biggest eigenvalues to be calculated
*
*   Outputs:
*       - double *lamda: calculated k eigenvalues will be stored here as a vector
*       - double *V1: calculated eigenvector will be returned by deferring
*/

    double lambda_old = 0;
    double num;
    double den;

    double* V0 = malloc(sizeof(double) * n);        // VO is nedded to find the corresponding eigenvector
    double* V1 = malloc(sizeof(double) * n);        // Corresponding eigenvector will be stored here


    for(int cont=0;cont<k;cont++){

        srand(time(NULL));
        // Initialize V0 with random numbers
        for(int i=0; i<n ; i++) V0[i] = rand()%100+1;


        for(int iter=1; iter<maxIter; iter++){


        for(int l=0;l<cont;l++){  //For que quita proyecciones

            double mag = dotProd(V0,FI[l],n);

            for(int m=0;m<n;m++){
                V0[m] -= mag*FI[l][m];
            }

        }
        ////////////////////////////////////////////////////
        normV(V0,n);            // Normalize V0 to avoid working with large values
        matXvec(A,V0,V1,n);     // We get V1 by multiplying A*V0

        double num = 0;
        double den = 0;

        // The value proposed for the eigenvalue is proposed
        for(int i=0; i<n; i++){
            num += V1[i]*V1[i];
            den += V1[i]*V0[i];
        }

        lamda[cont] = num/den;

        // lambda will be set as the eigenvalue only if the difference with the previous value is far too small,
        // this would mean that the value has converged
        if(fabs(lamda[cont]-lambda_old)<toler){
            normV(V1,n);                               // Before leaving this loop the eigenvector found is normalized
            for(int i=0;i<n;i++) FI[cont][i] = V1[i];
            break;                                    // break to leave this foor loop
        }
        else{                                   // In the case the difference between the current proposal to eigenvalue and the last one is significant
            lambda_old = lamda[cont];               // lambda_old will take the value of lambda to start next iteration
            for(int i=0;i<n;i++) V0[i] = V1[i];   // as VO will take values of V1
        }
        }
        ////////////////////////////////////////////////
    }

    free(V0);
    free(V1);
}
