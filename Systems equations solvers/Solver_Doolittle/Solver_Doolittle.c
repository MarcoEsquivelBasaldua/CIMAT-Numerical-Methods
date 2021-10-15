/*
*   Code name: Max_Index.c
*   Author: Marco Antonio Esquivel Basaldua
*   Date: 21/08/19
*
*   Code description:
*   This Code captures Matrix A,
*   captures the values of solutions to n linear equations in the matrix b,
*   If A is a square Matrix and b has sufficient values as solutions for the linear system:
*      calculates the values of x from the equation Ax = b, by the factorization of A = LU using Doolittle algorithm
*      get the determinant of the matrix A
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void Solver_Doolitle(double**,double**,int,int,int,int);    // This function factors A in LU and calculates the solutions to x given tha matrix b with the solutions to the n linear equations
int is_possible(double**,double**,int,int,int,int);        // Function to check if x could be found given A and b
double* solver_TSupAug(double**,int);                     // Function to solve x given an upper triangular augmented matrix
double* solver_TInfAug(double**,int);                    // Function to solve x given a lower triangular augmented matrix
double determinant(double**,int);                       // Function to get the determinant of a or triangular diagonal matrix

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

    // funtion to solve linear system
    Solver_Doolitle(A, b, n, m, k, l);

    free(A);
    free(b);
    return 0;
}

void Solver_Doolitle(double** A, double** b, int n, int m, int k, int l){
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
    double **L;                                                   // Lower triangular matrix factored
    L = malloc(n * sizeof *L);
    for (int i=0; i<n; i++) L[i] = malloc((m+1) * sizeof *L[i]);  // Matrrix L has an extra column so it could be transformed to an augmented matrix

    double **U;                                                   // Upper triangular matrix factored
    U = malloc(n * sizeof *U);
    for (int i=0; i<n; i++) U[i] = malloc((m+1) * sizeof *U[i]);  // Matrrix U has an extra column so it could be transformed to an augmented matrix

    double* x = malloc(sizeof(double) * n);           // Array where solution to x will be stored LUx = b
    double* y = malloc(sizeof(double) * n);          // Array where solution to y will be stored Ly = b
    int up = 0;                                     // Variable used to convert to upper triangular matrix
    int x_index[n];                                // In this vector indices of x order solution will be stored
    int nSwaps = 1;                               // Keep track if rows or columns are swapped
    double det_A;                                // Determinant of the matrix
    int flag;                                   // Flag will tell if it is possible to perform the calculations

    // First evaluates if this process can be done
    flag = is_possible(A,b,n,m,k,l);

    // Calculates and prints solutions to unknowns x and determinanat if possible
    if(flag == 1){

        // Initialize 1's and 0's in L and U
        for(int i=0;i<n;i++){
            for(int j=0;j<n;j++){
                if(i<j) L[i][j] = 0.0;
                else if(i>j) U[i][j] = 0.0;
                else L[i][j] = 1.0;
            }
        }


        int swap_occured = 0; // If a swap between rows is nedded this counter tells how many have accurred

        for(int j=0;j<n;j += 1-swap_occured){
            double sum = 0;
            swap_occured = 0;

            // First row of U and first column of L are calculated
            if(j==0){
                for(int j=0;j<n;j++) U[0][j] = A[0][j];
                for(int i=1;i<n;i++) L[i][0] = A[i][0]/U[0][0];
            }

            // Use of Doolittle Algorithm
            for(int i=1; i<=j;i++){
                for(int k=0;k<=i-1;k++){
                    sum += L[i][k]*U[k][j];
                }
                U[i][j] = A[i][j]-sum;
                sum = 0;
            }

            for(int i=j+1;i<n;i++){
                for(int k=0;k<=j-1;k++){
                    sum += L[i][k]*U[k][j];
                }
                L[i][j] = (A[i][j]-sum)/U[j][j];
                sum = 0;
            }
            ///////////////////////

            // Swap between rows if element on U[j][j]==0
            if(U[j][j] == 0){
                nSwaps *= -1;
                swap_occured++;

                // Swap of rows in matrix A
                double row_aux[n+1];
                for(int i=0;i<n;i++){
                    row_aux[i] = A[j+1][i];
                    A[j+1][i] = A[j][i];
                    A[j][i] = row_aux[i];
                }

                // Swap of rows in matrix U
                for(int i=0;i<n;i++){
                    row_aux[i] = U[j+1][i];
                    U[j+1][i] = U[j][i];
                    U[j][i] = row_aux[i];
                }

                // Swap of rows in matrix L
                for(int i=0;i<n;i++){
                    row_aux[i] = L[j+1][i];
                    L[j+1][i] = L[j][i];
                    L[j][i] = row_aux[i];
                }

                // Swap of rows in matrix b
                double b_aux= b[j+1][0];
                b[j+1][0] = b[j][0];
                b[j][0] = b_aux;

                // 0's and 1's are restored
                for(int i=0;i<n;i++){
                    for(int j=0;j<n;j++){
                        if(i<j) L[i][j] = 0.0;
                        else if(i>j) U[i][j] = 0.0;
                        else L[i][j] = 1.0;
                    }
                }
            }
        }


        // Print L and U
        printf("\nL is equal to:\n");
        for(int i=0;i<n;i++){
            for(int j=0;j<n;j++) printf("%.10lf ",L[i][j]);
            printf("\n");
        }

        printf("\nU is equal to:\n");
        for(int i=0;i<n;i++){
            for(int j=0;j<n;j++) printf("%.10lf ",U[i][j]);
            printf("\n");
        }



        // Calculates the determinant of A calculating the determinant of U
        det_A = nSwaps*determinant(U,n); //pow(-1.0,nSwaps)*

        // If the determinant could be calculated and it is different from 0 it means x could be calculated
        if(isnan(det_A) || det_A == 0){
            printf("Error solving x: check if matrix A is invertible or try rearranging its rows\n\n");

        }
        else
        {
            // L is set as a lower triangular augmented matrix with b as solutions
            for(int i=0;i<n;i++) L[i][m] = b[i][0];
            // y is solved from Ly = b
            y = solver_TInfAug(L,n);

            // U is set as an upper triangular augmented matrix with y as solutions
            for(int i=0;i<n;i++) U[i][m] = y[i];
            // x is solved from LUx = b
            x = solver_TSupAug(U,n);

            printf("\nThe values of the unknowns are:\n");
            for(int i=0;i<n;i++) printf("%lf\n",x[i]);
            printf("\nThe determinant of the matrix A is:\n%lf\n\n",det_A);
        }
    }

    free(L);
    free(U);
    free(x);
    free(y);
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

double* solver_TSupAug  (double** A, int n){
    /*
*   Inputs:
*       - double** A: pointer to augmented matrix
*       - int n: mumber of rows of matrix A
*
*   Output:
*       - double* x: x solutions to augmented upper matrix A
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

double* solver_TInfAug  (double** A, int n){
    /*
*   Inputs:
*       - lomg double** A: pointer to augmented matrix
*       - int n: mumber of rows of matrix A
*
*   Output:
*       - long double* x: x solutions to augmented lower matrix A
*/
    double* x = malloc(sizeof(double) * n);
    double sum = 0;

    x[0] = A[0][n]/A[0][0];

    for(int i=1;i<n;i++){
            sum = 0;
            for(int j=0;j<i;j++) sum += A[i][j]*x[j];
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
