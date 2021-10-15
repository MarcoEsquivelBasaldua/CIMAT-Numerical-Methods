/*
*   Code name: solver_Gauss.c
*   Author: Marco Antonio Esquivel Basaldua
*   Date: 21/08/19
*
*   Code description:
*   This Code captures Matrix A,
*   captures the values of solutions to n linear equations in the vector b,
*   If A is a square Matrix and b has sufficient values as solutions for the linear system:
*      calculates the values of x from the equation Ax = b
*      get the determinant of the matrix A
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void solver_Gauss(double**,double**,int,int,int,int);    // Function to solve x given a square matrix ussing Gauss elimination
int is_possible(double**,double**,int,int,int,int);    // Function to check if x could be found given A and b
double determinant(double**,int);                     // Function to get the determinant of a or triangular diagonal matrix
double** toUpTriMatrix(double**,int );              // Transforms a given augmented matrix to upper triangular
double* solver_TSpuAug(double**,int);             // Function to solve x given an upper triangular augmented matrix


int main(int argc , char* argv[] ){
    int n,m;        // Size of matrix A
    int k,l;        // Size of matrix b (solutions to linear equations)
    double **A;     // Matrix A
    double **b;     // Values of solutions to the n linear equations
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
    for (i=0; i<n; i++) A[i] = malloc((m+1) * sizeof *A[i]);  // Matrrix A has an extra column so it could be transformed to an augmented matrix

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
    solver_Gauss(A, b, n, m, k, l);

    free(A);
    free(b);
    return 0;
}

void solver_Gauss(double** A, double** b, int n, int m, int k, int l){
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
        // Matrices A and b are stored as only one augmented matrix
        for(i=0; i<n; i++) A[i][m] = b[i][0];

        // New Augmented matrix A is set to upper triangular matrix
        A = toUpTriMatrix(A,n);

        // Determinanat of A is calculated
        det_A = determinant(A,n);

        // If the determinant could be calculated and it is different from 0 it means x could be calculated
        if(isnan(det_A) || det_A == 0){
            printf("Error solving for x: check if matrix A is invertible or try rearranging its rows\nUsing Solver_GaussPivot is recomended.\n\n\n");

        }
        else
        {
            // Solve for x using a new version of solver_TSup wich only requires the upper triangular aumented matrix and number of rows
            x = solver_TSpuAug(A,n);

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

double** toUpTriMatrix(double** A,int n){
/*
*   Inputs:
*       - double** A: pointer to augmented matrix
*       - int n: mumber of rows of matrix A
*
*   Outputs:
*       - double** A: Upper triangular matrix equivalent for input A
*/
    int upp = 0;    // This variable kepps track on how many elements in the column have been set to zero by iteration
    for(int q=0;q<=n;q++){
        upp = 0;
        for(int p=q+1;p<n;p++){
            for(int r=q+1;r<=n;r++){
                A[p][r] = A[p][r]-(A[p-1-upp][r]/A[q][q])*A[p][q];
            }
        A[p][q]=0;
        upp++;
        }
    }

    return A;
}

double* solver_TSpuAug  (double** A, int n){
/*
*   Inputs:
*       - double** A: pointer to augmented matrix
*       - int n: mumber of rows of matrix A
*
*   Output:
*       - double* x: x solutions to augmented matrix A
*/
    double* x = malloc(sizeof(double) * n);
    double sum = 0;

    x[n-1] = A[n-1][n]/A[n-1][n-1];

    for(int i=n-2;i>=0;i--){
            sum = 0;
            for(int j=n-1;j>i;j--) sum += A[i][j]*x[j];
            x[i] = (A[i][n] - sum)/A[i][i];
    }

    return x;
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
