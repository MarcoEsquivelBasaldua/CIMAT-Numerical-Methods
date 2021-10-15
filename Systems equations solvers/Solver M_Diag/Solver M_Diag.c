/*
*   Code name: solver_MDiag.c
*   Author: Marco Antonio Esquivel Basaldua
*   Date: 25/08/19
*
*   Code description:
*   This Code captures Matrix A,
*   captures the values of solutions to n linear equations in the matrix b,
*   If A is a square diagonal Matrix and b has sufficient values as solutions for the linear system:
*      calculates the values of x from the equation Ax = b
*      get the determinant of the matrix A
*/
#include <stdio.h>
#include <stdlib.h>

void solver_MDiag(double**,double**,int,int,int,int);    // Function to solve x given a diagonal matrix A and solutions to the linear equations b
int is_possible(double**,double**,int,int,int,int);    // Function to check if x could be found given A and b
double determinant(double**,int);                     // Function to get the determinant of a diagonal matrix A


int main( int argc , char* argv[] ){
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

    // Funtion to solve x in a linear system
    solver_MDiag(A, b, n, m, k, l);


    free(A);
    free(b);
    return 0;
}

void solver_MDiag(double** A, double** b, int n, int m, int k, int l){
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
    double* x = malloc(sizeof(double) * n);     // Array where solutions to x will be stored
    double det_A;                              // Variable to store the determinant of A
    int flag;                                 // If flag is equal to 1, solutions to x are solved

    // First, evaluates if this process can be done
    flag = is_possible(A,b,n,m,k,l);

    // Calculates and prints solutions to unknowns x and determinant if possible
    if(flag == 1){
        printf("The values of the unknowns are:\n");
        for(int i=0;i<n;i++){
            x[i] = b[i][0]/A[i][i];
            printf("%lf\n",x[i]);
        }

        det_A = determinant(A,n);
        printf("\nThe determinant of the matrix A is:\n%lf\n\n",det_A);
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
    int flag = 1;       // flag = 1 if solutions to x can be found

    // Check if A is a square matrix
    if(n != m){
        flag = 0;
        printf("\nERROR: A is not a square matrix\n");
    }

    // Check if A is a diagonal matrix
    if(flag == 1){
        for(int i=0; i<n; i++){
            for(int j=0; j<m; j++){
                if(i == j){
                    if(A[i][j] == 0){
                        flag = 0;
                        break;
                    }
                }
                else{
                    if(A[i][j] != 0){
                        flag = 0;
                        break;
                    }
                }
            }
            if(flag == 0) break;
        }
        if(flag == 0) printf("\nERROR: A is not a diagonal matrix\n");
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
