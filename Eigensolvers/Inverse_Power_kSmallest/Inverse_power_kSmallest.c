/*
*   Code name: Inverse_power_kSmallest.c
*   Author: Marco Antonio Esquivel Basaldua
*   Date: 14/09/19
*
*   Code description:
*   This Code captures Matrix A and finds its k smallest eigenvalues and their
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

void Inverse_Power(double**,double*,double**,int,double,int,int);       // Inverse Power method to find k smallest eigenvalues and eigenvectors assotiated
void LU_Fact(double**,double**,double**,int);                          // This function factors a given squared matrix A as A = LU
void normV(double*,int);                                              // Finds the euclidean norm from a vector
void solver_TSupAug(double**,int,double*);                           // Function to solve x given an upper triangular augmented matrix
void solver_TInfAug(double**,int,double*);                          // Function to solve x given a lower triangular augmented matrix
double dotProd(double*,double*,int);                               // Calculates the dot product given two vectors
void Resuelve(double**,double**,double*,double*,int);             // Finds the values of the unknowns given a linear sistem Ax=b wich A matrix is factored in L and U

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
    FI = calloc(K, sizeof *FI);
    for (int i=0; i<K; i++) FI[i] = malloc(n * sizeof *FI[i]);


    Inverse_Power(A,lamda,FI,n,toler,maxIter,K);

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

void Resuelve(double **L,double **U,double *V0,double *V1,int n){
/*
*   Intputs:
*       - double **L: L factor from a previous A matrix
*       - double **U: U factor from a previous A matrix
*       - double *b: solutions to the n linear equations in Ax = b
*       - double *x: vector with values of the unknowns, this vector is passed by deferring
*       - int n: number of columns in L, U, b and x
*
*   Outputs:
*       - double *x: vector with values of the unknowns, this vector is passed by deferring
*/
    double* y = malloc(sizeof(double) * n);           // Array where solution to x will be stored LUx = b

    for(int i=0;i<n;i++) L[i][n] = V0[i];
    // y is solved from Ly = b
    solver_TInfAug(L,n,y);

    // U is set as an upper triangular augmented matrix with y as solutions
    for(int i=0;i<n;i++) U[i][n] = y[i];
    // x is solved from LUx = b
    solver_TSupAug(U,n,V1);

    free(y);
}

void LU_Fact(double **A,double **L,double **U, int n){
/*
*   Inputs:
*       - double **A: given sqared matrix to be factored as A = LU
*       - double **L: L factor from A matrix, this matrix will be returnen by deferring
*       - double **U: U factor from A matrix, this matrix will be returnen by deferring
*       - int n: dimension of square matrix A
*
*   Outputs:
*       - double **L: L factor from A matrix, this matrix will be returnen by deferring
*       - double **U: U factor from A matrix, this matrix will be returnen by deferring
*/
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
        long double sum = 0;

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
        }
        else{
            swap_occured = 0;
            swap=0;
            j++;
        }
    }
}

void solver_TSupAug  (double** A, int n, double* x){
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

void solver_TInfAug  (double** A, int n, double* x){
    /*
*   Inputs:
*       - lomg double** A: pointer to augmented matrix
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

void Inverse_Power(double** A,double* lamda,double** FI,int n,double toler,int maxIter,int k){
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
    double lambda_old = 100000;
    double num;
    double den;

    double* V0 = malloc(sizeof(double) * n);        // VO is nedded to find the corresponding eigenvector
    double* V1 = malloc(sizeof(double) * n);        // Corresponding eigenvector will be stored here

    double **L;                                                   // Lower triangular matrix factored
    L = calloc(n, sizeof *L);
    for (int i=0; i<n; i++) L[i] = malloc((n+1) * sizeof *L[i]);  // Matrrix L has an extra column so it could be transformed to an augmented matrix

    double **U;                                                   // Upper triangular matrix factored
    U = calloc(n , sizeof *U);
    for (int i=0; i<n; i++) U[i] = malloc((n+1) * sizeof *U[i]);  // Matrrix U has an extra column so it could be transformed to an augmented matrix

    LU_Fact(A,L,U,n);

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

            Resuelve(L,U,V0,V1,n);

            double num = 0;
            double den = 0;

            // The value proposed for the eigenvalue is proposed
            for(int i=0; i<n; i++){
                num += V1[i]*V0[i];
                den += V1[i]*V1[i];
            }

            lamda[cont] = num/den;

            // lambda will be set as the eigenvalue only if the difference with the previous value is far too small,
            // this would mean that the value has converged
            if(fabs(lamda[cont]-lambda_old)<toler){
                normV(V0,n);                               // Before leaving this loop the eigenvector found is normalized
                for(int i=0;i<n;i++) FI[cont][i] = V0[i];
                break;                                    // break to leave this foor loop
            }
            else{                                   // In the case the difference between the current proposal to eigenvalue and the last one is significant
                lambda_old = lamda[cont];               // lambda_old will take the value of lambda to start next iteration
                for(int i=0;i<n;i++) V0[i] = V1[i];   // as VO will take values of V1
            }
        }
        ////////////////////////////////////////////////
    }

    for(int i=0;i<n;i++) free(L[i]);
    free(L);
    for(int i=0;i<n;i++) free(U[i]);
    free(U);

    free(V0);
    free(V1);
}
