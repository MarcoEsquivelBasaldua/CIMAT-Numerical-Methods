/*
*   Code name: solver_TSup.c
*   Author: Marco Antonio Esquivel Basaldua
*   Date: 25/08/19
*
*   Code description:
*   This Code captures Matrix A,
*   captures the values of solutions to n linear equations in the matrix b,
*   If A is a square upper triangular Matrix and b has sufficient values as solutions for the linear system:
*      calculates the values of x from the equation Ax = b
*      get the determinant of the matrix A
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void solver_TSup(double**,double**,int,int,int,int);    // Function to solve x given a square upper triangular matrix
int is_possible(double**,double**,int,int,int,int);    // Function to check if x could be found given A and b
double determinant(double**,int);                     // Function to get the determinant of a diagonal or triangular matrix


int main(int argc , char* argv[] ){
    int n,m;        // Size of matrix A
    int k,l;        // Size of matrix b (solutions to linear equations)
    double **A;     // Matrix A
    double **b;     // Values of solutions to the n linear equations
    int i,j;        // Integers used in all for cycles in main

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

    for(i=0; i<n; i++){
		for(j=0; j<m; j++){
			fscanf(fin, "%lf", &A[i][j]);
			printf("%lf ", A[i][j]);
		}
		printf("\n");
	}
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

    for(i=0; i<k; i++){
		for(j=0; j<l; j++){
			fscanf(fin, "%lf", &b[i][j]);
			printf("%lf ", b[i][j]);
		}
		printf("\n");
	}
    printf("\n");
	fclose( fin );

    // Funtion to solve linear system
    solver_TSup(A, b, n, m, k, l);


    free(A);
    free(b);
    return 0;
}

void solver_TSup(double** A, double** b, int n, int m, int k, int l){
/*
*   Inputs:
*       - double **A: Pointer to the matrix A
*       - double **b: Pointer to the matrix b
*       - int n: number of rows in A
*       - int m: number of columns in A
*       - int k: number of rows in b
*       - int l: number of columns in b
*
*   Outputs:
*       - None
*/
    double* x = malloc(sizeof(double) * n);     // Vector where solutions to x will be stored
    double det_A;                              // Variable to store the determinant of A
    int flag;                                 // If flag is equal to 1, solutions to x are solved
    double sum = 0;                          // This variable sum will be used acording to the algorithm sol dolve for x
    int i,j;                                // Indices used in for cycles

    // First evaluates if this process can be done
    flag = is_possible(A,b,n,m,k,l);

    // Calculates and prints solutions to unknowns x and determinanat if possible
    if(flag == 1){
        // Determinanat of A is calculated
        det_A = determinant(A,n);

        // If the determinant could be calculated and it is different from 0 it means x could be calculated
        if(isnan(det_A) || det_A == 0){
            printf("Error solving for x: check if matrix A is invertible or try rearranging its rows\n\n");
        }
        else
        {
            x[n-1] = b[n-1][0]/A[n-1][n-1];

            for(i=n-2;i>=0;i--){
                sum = 0;
                for(j=n-1;j>i;j--) sum += A[i][j]*x[j];
                x[i] = (b[i][0] - sum)/A[i][i];
            }

            printf("\nThe values of the unknowns are:\n");
            for(i=0;i<n;i++) printf("%lf\n",x[i]);
            printf("\nThe determinant of the matrix A is:\n%lf\n\n",det_A);
        }
    }

    free(x);
}

int is_possible(double** A, double** b, int n, int m, int k, int l){
/*
*   Inputs:
*       - double **A: Pointer to the matrix A
*       - double **b: Pointer to the matrix b
*       - int n: number of rows in A
*       - int m: number of columns in A
*       - int k: number of rows in b
*       - int l: number of columns in b
*
*   Outputs:
*       - None
*/
    int flag = 1;       // flag = 1 if solutions to x can be found, flag = 0 otherwise

    // Check if A is a square matrix
    if(n != m){
        flag = 0;
        printf("\nERROR: A is not a square matrix\n");
    }

    // Check if A is an upper triangular matrix
    if(flag == 1){
        for(int i=0; i<n; i++){
            for(int j=0; j<m; j++){
                if(i == j){
                    if(A[i][j] == 0){
                        flag = 0;
                        break;
                    }
                }
                else if(i > j){
                    if(A[i][j] != 0){
                        flag = 0;
                        break;
                    }
                }
            }
            if(flag == 0) break;
        }
        if(flag == 0) printf("\nERROR: A is not an upper triangular matrix\n");
    }

    // Check if b is a vector
    if(flag == 1 && l != 1 )
    {
        flag = 0;
        printf("\nERROR: b must be a column vector\n");
    }

    // Check if b has enough values as rows in A
    if(flag == 1 && n != k )
    {
        flag = 0;
        printf("\nERROR: b must provide as many values as rows in matrix A\n");
    }

    return flag;
}

double determinant(double **A, int n){
/*
*   Inputs :
*       - double **A: pointer to square matrix A
*       - int n: number of rows (or columns) of matrix A
*
*   Outputs:
*       - double det_A: determinant of matrix A
*/
    double det_A = 1.0;
    for(int i=0;i<n;i++) det_A *= A[i][i];

    return det_A;
}
