/*
*   Code name: Rayleigh.c
*   Author: Marco Antonio Esquivel Basaldua
*   Date: 21/09/19
*
*   Code description:
*   This Code captures Matrix A and finds one of its eigenvalues and its eigenvector
*   assotiated using the Rayleigh method
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

double Rayleigh(double**,double*,double*,int,double,int);    // Rayleigh method to find one of its eigenvalues and eigenvector assotiated
double Ro_quotient(double*,double**,int);                   // Calculates the Rayleigh quotient
void Solver_Doolitle(double**,double*,double*,int);        // This function factors A in LU and calculates the solutions to x given tha matrix b with the solutions to the n linear equations
void solver_TSupAug(double**,int,double*);                // Function to solve x given an upper triangular augmented matrix Ux = b
void solver_TInfAug(double**,int,double*);               // Function to solve x given a lower triangular augmented matrix Lx = b
double dotProd(double*,double*,int);                    // Calculates the dot product given two vectors
void matXvec(double**,double*,double*,int);            // Function to multiply a matrix times a vector
void V_normalize(double*,int);                        // normalizes a vector
double normV(double*,int);                           // Finds the euclidean norm of a vector

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


    double* eV = malloc(sizeof(double) * n);            // Corresponding eigenvector will be stored here
    double* eV_approx = malloc(sizeof(double) * n);    // An approximation to the eigenvector

    double lamda = 0;                         // Found eigenvalue will be stored here
    double toler = 1e-12;                    // A tolerance will determine if the code has found the solutions
    int maxIter = 10000;                    // This value is the maximum of iterations the code will run. This avoids the code to execute forever

    srand(time(NULL));
    // Initialize eV_approx with random numbers
    for(i=0; i<n ; i++) eV_approx[i] = rand()%100+1;

    // Time excecution starts
    clock_t t_ini, t_fin;
    double secs;
    t_ini = clock();

    lamda = Rayleigh(A,eV,eV_approx,n,toler,maxIter);
    printf("\n");

    // Time excecution is displayed
    t_fin = clock();
    secs = (double)(t_fin - t_ini) / CLOCKS_PER_SEC;
    printf("%.16g miliseconds taken in the algorithm\n", secs * 1000.0);

    // Memory is freed
    for(int i=0;i<n;i++) free(A[i]);
    free(A);
    free(eV);
    free(eV_approx);

    return 0;
}

double Rayleigh(double **A,double *eVnew, double *eVold,int n,double toler,int maxIter){
/*
*   Inputs:
*       - double **A: square matrix to calculate one of its eigenvalues and its assotiated eigenvector
*       - double *eVnew: calculated eigenvector will be stored here
*       - double *eVold: initialized eigenvector approximation
*       - int n: number of rows in matrix A and vectors eVnew and eVold
*       - double toler: tolerance near to cero to determine if the code has found the solutions
*       - int max_iter: This value is the maximum of iterations the code will run. This avoids the code to execute forever
*
*   Outputs:
*       - double lamda: found eigenvalue of matrix A
*       - double *eVnew: calculated eigenvector will be returned by deferring
*/

    double **Acopy;     // A copy of matrix A is made this will avoid to change the original matrix

    Acopy = malloc(n * sizeof *Acopy);
    for (int i=0; i<n; i++) Acopy[i] = malloc(n * sizeof *Acopy[i]);

    for(int i=0; i<n; i++) for(int j=0; j<n; j++) Acopy[i][j] = A[i][j];
    /////

    double* Ay = malloc(sizeof(double) * n);    // This vector will help to test the termination condition of the

    V_normalize(eVold,n);

    double lambda_old = Ro_quotient(eVold,A,n);
    double lambda;

    // Rayleigh algorithm begins
    for(int iter=1; iter<maxIter; iter++){

        for(int i=0; i<n; i++) Acopy[i][i] = A[i][i]-lambda_old;

        Solver_Doolitle(Acopy,eVold,eVnew,n);

        V_normalize(eVnew,n);

        lambda = Ro_quotient(eVnew,A,n);

        matXvec(A,eVnew,Ay,n);

        for(int i=0; i<n; i++) Ay[i] = Ay[i]-lambda*eVnew[i];

        double norm = normV(Ay,n);

        if(norm < toler){
            printf("Iterations taken: %d\n", iter);     // This will help know how many iterations it took to get to the results
            break;                                    // break to leave this foor loop
        }
        else{                                   // In the case the difference between the current proposal to eigenvalue and the last one is significant
            lambda_old = lambda;               // lambda_old will take the value of lambda to start next iteration
            for(int i=0;i<n;i++) eVold[i] = eVnew[i];   // as eVold will take values of eVnew
        }
    }

    // Results are displayed
    printf("Found Eigenvalue:\n %g\n", lambda);

    printf("Eigenvector assotiated to found eigenvalue is:\n");
    for(int i=0; i<n; i++) printf("%.10g\n",eVnew[i]);



    for(int i=0;i<n;i++) free(Acopy[i]);
    free(Acopy);
    free(Ay);
    return lambda;
}

// Quotient of Rayleigh
double Ro_quotient(double *u,double **A,int n){
/*
*   Inputs:
*       - double *u
*       - double **A
*       - int n: number of rows in matrix A and vector u
*
*   Outputs:
*       - double : Rayleigh quotient = (u'Au)/(u'u)
*/
    double* Au = malloc(sizeof(double) * n);        // A*u multiplication
    double utAu;
    double utu;

    matXvec(A,u,Au,n);

    utAu = dotProd(u,Au,n);
    utu = dotProd(u,u,n);

    free(Au);

    return utAu/utu;
}

void Solver_Doolitle(double** A, double* b, double* x, int n){
/*
*   Inputs:
*       - double **A: Pointer to the matrix A
*       - double *b: Pointer to the vector b
*       - double *x: Pointer to the vector x (solutions to unknowns)
*       - int n: number of rows in A, b and x
*
*   Outputs:
*       - double *x: Pointer to the vector x (solutions to unknowns) by deferring
*/
    double **L;                                                   // Lower triangular matrix factored
    L = malloc(n * sizeof *L);
    for (int i=0; i<n; i++) L[i] = malloc((n+1) * sizeof *L[i]);  // Matrrix L has an extra column so it could be transformed to an augmented matrix

    double **U;                                                   // Upper triangular matrix factored
    U = malloc(n * sizeof *U);
    for (int i=0; i<n; i++) U[i] = malloc((n+1) * sizeof *U[i]);  // Matrrix U has an extra column so it could be transformed to an augmented matrix

    //double* x = malloc(sizeof(double) * n);           // Array where solution to x will be stored LUx = b
    double* y = malloc(sizeof(double) * n);          // Array where solution to y will be stored Ly = b
    int up = 0;                                     // Variable used to convert to upper triangular matrix
    int x_index[n];                                // In this vector indices of x order solution will be stored




    // Initialize 1's and 0's in L and U
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            if(i<j) L[i][j] = 0.0;
            else if(i>j) U[i][j] = 0.0;
            else L[i][j] = 1.0;
        }
    }


    int swap_occured = 0; // If a swap between rows is nedded this counter tells how many have accurred
    int swap = 0;
    int j=1;

    // First row of U and first column of L are calculated
    for(int j=0;j<n;j++) U[0][j] = A[0][j];
    for(int i=1;i<n;i++) L[i][0] = A[i][0]/U[0][0];

    while(j<n){//for(int j=0;j<n;j += 1-swap){
        double sum = 0;

        for(int i=1; i<=j; i++){
            for(int k=0;k<=i-1;k++){
                sum=sum+(L[i][k]*U[k][j]);
            }
            U[i][j]=A[i][j]-sum;
            sum=0;
        }

        for(int i=j+1;i<n;i++){
            for(int k=0;k<=j-1;k++){
                sum=sum+(L[i][k]*U[k][j]);
            }
            L[i][j]=(A[i][j]-sum)/U[j][j];
            sum=0;
        }

        // Swap between rows if element on U[j][j]==0
        if(U[j][j] == 0){
            swap_occured++;
            swap=1;

            // Swap of rows in matrix A
            double row_aux[n+1];
            for(int i=0;i<n;i++){
                row_aux[i] = A[j+swap_occured][i];
                A[j+swap_occured][i] = A[j][i];
                A[j][i] = row_aux[i];
            }

            // Swap of rows in matrix b
            double b_aux= b[j+swap_occured];
            b[j+swap_occured] = b[j];
            b[j] = b_aux;

        }
        else{
            swap_occured = 0;
            swap=0;
            j++;
        }
    }



    // L is set as a lower triangular augmented matrix with b as solutions
    for(int i=0;i<n;i++) L[i][n] = b[i];
    // y is solved from Ly = b
    solver_TInfAug(L,n,y);

    // U is set as an upper triangular augmented matrix with y as solutions
    for(int i=0;i<n;i++) U[i][n] = y[i];
    // x is solved from LUx = b
    solver_TSupAug(U,n,x);


    //printf("\nThe values of the unknowns are:\n");
    //for(int i=0;i<n;i++) printf("%lf\n",x[i]);


    for(int i=0;i<n;i++) free(L[i]); // Se libera la memoria solicitada
    free(L);
    for(int i=0;i<n;i++) free(U[i]); // Se libera la memoria solicitada
    free(U);
    free(y);
}


void solver_TSupAug (double** A, int n, double* x){
    /*
*   Inputs:
*       - double** A: pointer to augmented matrix
*       - int n: mumber of rows of matrix A
*
*   Output:
*       - double* x: x solutions to augmented upper matrix A
*/

    double sum = 0;

    x[n-1] = A[n-1][n]/A[n-1][n-1];

    for(int i=n-2;i>=0;i--){
            sum = 0;
            for(int j=n-1;j>i;j--) sum += A[i][j]*x[j];
            x[i] = (A[i][n] - sum)/A[i][i];
    }

}

void solver_TInfAug (double** A, int n, double* x){
    /*
*   Inputs:
*       - double** A: pointer to augmented matrix
*       - int n: mumber of rows of matrix A
*
*   Output:
*       - long double* x: x solutions to augmented lower matrix A
*/

    double sum = 0;

    x[0] = A[0][n]/A[0][0];

    for(int i=1;i<n;i++){
            sum = 0;
            for(int j=0;j<i;j++) sum += A[i][j]*x[j];
            x[i] = (A[i][n] - sum)/A[i][i];
    }

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

// This function normalizes a vector
void V_normalize(double *V,int n){
/*
*   Inputs:
*       - double *V: vector to be normalized, its corresponding vector normalized will be passed in this very same location
*       - int n: number of rows in vector V
*
*   Outputs:
*       - None
*/
    double norm = normV(V,n);

    for(int i=0; i<n; i++) V[i] = V[i]/norm;
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
