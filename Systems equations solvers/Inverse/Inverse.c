/*
*   Code name: Inverse.c
*   Author: Marco Antonio Esquivel Basaldua
*   Date: 21/08/19
*
*   Code description:
*   This Code captures Matrix A,
*   If A is a square Matrix and its determinant is different from 0
*      calculates the Inverse matrix of A
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double** Inverse(double **,int, int);                      // Function to get inverse Matrix of A
double* solver_TSupAug(double**,int);                     // Function to solve x given an upper triangular augmented matrix
double* solver_TInfAug(double**,int);                    // Function to solve x given a lower triangular augmented matrix
double determinant(double**,int);                       // Function to get the determinant of a or triangular diagonal matrix
void AxB(double**,int,int,double**,int,int);

int main(int argc , char* argv[] ){
    int n,m;        // Size of matrix A
    double **A;     // Matrix A
    int i,j;        // Integers used in all for cycles

    double **A_inv;     // Matrix A_inv where inverse of A will be stored
    A_inv = malloc(n * sizeof *A_inv);
    for (int i=0; i<n; i++) A_inv[i] = malloc(m * sizeof *A_inv[i]);

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

    // Function to calcul and return the inverse of A
    A_inv = Inverse(A,n,m);

    // A and A_inv are multiplied to demonstrate A*A_inv = I
    AxB(A,n,m,A_inv,n,m);

    printf("\n\n");

    free(A);
    free(A_inv);
    return 0;
}

double** Inverse(double **A, int n, int m){
/*
*   Inputs:
*       - double **A: Pointer to the matrix A
*       - int n: number of rows in A
*       - int m: number of columns in A
*
*   Outputs:
*       - double **A_inv: Inverse matrix of A
*/
    double **A_aux;     // Matrix A_aux is used as a copy of A, so A will reamin unchanged
    A_aux = malloc(n * sizeof *A_aux);
    for (int i=0; i<n; i++) A_aux[i] = malloc(m * sizeof *A_aux[i]);
    for(int i=0; i<n;i++)  for(int j=0; j<n; j++) A_aux[i][j] = A[i][j];

    double **A_inv;     // Matrix A_inv where inverse of A will be stored
    A_inv = malloc(n * sizeof *A_inv);
    for (int i=0; i<n; i++) A_inv[i] = malloc(m * sizeof *A_inv[i]);

    double **L;
    L = malloc(n * sizeof *L);                                          // Lower triangular matrix factored
    for (int i=0; i<n; i++) L[i] = malloc((n+1) * sizeof *L[i]);       // Matrrix L has an extra column so it could be transformed to an augmented matrix

    double **U;
    U = malloc(n * sizeof *U);                                          // Upper triangular matrix factored
    for (int i=0; i<n; i++) U[i] = malloc((n+1) * sizeof *U[i]);       // Matrrix U has an extra column so it could be transformed to an augmented matrix

    int flag;                                   // Flag will tell if it is possible to perform the calculations
    int nSwaps = 1;                            // Keep track if rows or columns are swapped
    double det_A;                             // Determinant of the matrix

    // A_inv is first set to the identity matrix
    for(int i=0; i<n;i++){
        for(int j=0; j<n; j++){
            if(i == j) A_inv[i][j] = 1.0;
            else A_inv[i][j] = 0.0;
        }
    }

    // First evaluates if this process can be done, A has to be a squared matrix
    if(n == m){
        // Get LU factorization from A
        // Initialize 1's and 0's in L and U
        for(int i=0;i<n;i++){
            for(int j=0;j<n;j++){
                if(i<j) L[i][j] = 0.0;
                else if(i>j) U[i][j] = 0.0;
                else L[i][j] = 1.0;
            }
        }

        int swap_occured = 0;   // If a swap between rows is nedded this counter tells how many have accurred

        for(int j=0;j<n;j += 1+swap_occured){
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
                U[i][j] = A_aux[i][j]-sum;
                sum = 0;
            }

            for(int i=j+1;i<n;i++){
                for(int k=0;k<=j-1;k++){
                    sum += L[i][k]*U[k][j];
                }
                L[i][j] = (A_aux[i][j]-sum)/U[j][j];
                sum = 0;
            }
            //////////////////////////////

            // Swap between rows if element on U[j][j]==0
            if(U[j][j] == 0){
                nSwaps *= -1;
                swap_occured--;

                // Swap of rows in matrix A
                double row_aux[n+1];
                for(int i=0;i<n;i++){
                    row_aux[i] = A_aux[j+1][i];
                    A_aux[j+1][i] = A_aux[j][i];
                    A_aux[j][i] = row_aux[i];
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

                // Swap of rows in A_inv
                for(int i=0;i<n;i++){
                    row_aux[i] = A_inv[j+1][i];
                    A_inv[j+1][i] = A_inv[j][i];
                    A_inv[j][i] = row_aux[i];
                }

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

        printf("\n\nThe determinanat of the matrix A is:\n%lf\n",det_A);

        // If the determinant could be calculated and it is different from 0 it means x could be calculated
        if(isnan(det_A) || det_A == 0){
            printf("\n\nError solving x: check if matrix A is invertible or try rearranging its rows\n\n");

        }
        else
        {
            //For each column in A_inv = I
            double* x = malloc(sizeof(double) * n);           // Array where solution to x will be stored LUx = b
            double* y = malloc(sizeof(double) * n);           // Array where solution to y will be stored Ly = b
            for(int j=0;j<n;j++){

                for(int i = 0; i<n; i++) L[i][n] = A_inv[i][j];

                y = solver_TInfAug(L,n);

                // Take y as solution to augmented matrix Ux=y and solve for x
                for(int i = 0; i<n; i++) U[i][n] = y[i];

                x = solver_TSupAug(U,n);

                // Replace x in column j of A_inverse
                for(int i =0;i<n;i++) A_inv[i][j] = x[i];
            }

            // Shows A_inv
            printf("\n\nInverse of A is:\n");
            for(int i=0; i<n;i++){
                for(int j=0; j<n; j++){
                    printf("%lf ",A_inv[i][j]);
                }
                printf("\n");
            }

        }
    }

    free(L);
    free(U);
    free(A_aux);

    return A_inv;
}

double* solver_TSupAug  (double** A, int n){
/*
*   Inputs:
*       - lomg double** A: pointer to augmented matrix
*       - int n: mumber of rows of matrix A
*
*   Output:
*       - long double* x: x solutions to augmented upper matrix A
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

void AxB(double** A, int n, int m, double** B, int p, int q){
/*
*   Inputs:
*       - double** A: pointer to matrix A
*       - int n: row dimension of A
*       - int m: cols dimension of A
*       - double** B: pointer to matrix B
*       - int p: row dimension of B
*       - int q: cols dimension of B
*
*   Outputs:
*       - double AxB: multiplication A*B if possible
*/

    if(m == p){
        double** AxB;
        AxB = malloc(n * sizeof *AxB);
        for (int i=0; i<n; i++) AxB[i] = malloc(q * sizeof *AxB[i]);

        printf("\nAxB=\n");
        for(int i=0; i<n; i++){
            for(int j=0; j<q; j++){
                AxB[i][j] = 0;
                for(int k=0; k<m; k++) AxB[i][j] += A[i][k] * B[k][j];

                printf("%lf ",AxB[i][j]);
            }
            printf("\n");
        }
    }
    else{
        printf("\nError: A and B are NOT compatible\n Multiplication can not be done\n");
    }
}
