/*
*   Code name: solver_GaussPivot.c
*   Author: Marco Antonio Esquivel Basaldua
*   Date: 21/08/19
*
*   Code description:
*   This Code captures Matrix A,
*   captures the values of solutions to n linear equations in the matricx b,
*   If A is a square Matrix and b has sufficient values as solutions for the linear system:
*      calculates the values of x from the equation Ax = b
*      get the determinant of the matrix A
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void solver_GaussPivot(double**,double**,int,int,int,int);    // Function to solve x given a square matrix ussing Gauss elimination with pivots
double* solver_TSupAug(double**,int);                       // Function to solve x given an upper triangular augmented matrix
double determinant(double**,int);                          // Function to get the determinant of a diagonal matrix
int is_possible(double**,double**,int,int,int,int);       // Function to check if x could be found given A and b


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

    // funtion to solve linear system
    solver_GaussPivot(A, b, n, m, k, l);

    free(A);
    free(b);
    return 0;
}

void solver_GaussPivot(double** A, double** b, int n, int m, int k, int l){
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
    double* x = malloc(sizeof(double) * n);          // Array where solution to x will be stored
    double* x_sort = malloc(sizeof(double) * n);    // Array where solutions to x, once rearranged will be stored
    double max = 0, aux = 0;                       // Variables to store maximum value found in matrix
    int ind_i,ind_j;                              // Indices where max value is located
    int up = 0;                                  // Variable used to convert to upper triangular matrix
    int x_index[n];                             // In this vector indices of x order solution will be stored
    int nSwaps = 1;                            // Keep track if rows or columns are swapped
    double det_A;                             // Determinant of the matrix
    int flag;                                // Flag will tell if it is possible to perform the calculations

    // First evaluates if this process can be done
    flag = is_possible(A,b,n,m,k,l);

    // Calculates and prints solutions to unknowns x and determinanat if possible
    if(flag == 1){

        // Matrices A and b are stored as only one augmented matrix
        for(int i=0; i<n; i++) A[i][m] = b[i][0];

        // Initialize order in indices of x solutions
        for(int i =0; i<n; i++) x_index[i] = i;
        /////////////////////////////////////////

        //Convert to Up triangular matrix using Pivot
        for(int q=0;q<n;q++){
            up = 0;
            max = 0;

            // Find maximum value in coefficient matrix
            for(int i=q;i<n;i++){
                for(int j=q;j<n;j++){
                    if(A[i][j] < 0) aux = -A[i][j];
                    else aux = A[i][j];

                 if(aux > max){
                        max = aux;
                        ind_i = i;
                        ind_j = j;
                    }
                }
            }
            /////////////////////////////////////////////

            // Swap rows
            double row_aux[n+1];
            if(q!=ind_i){
                nSwaps *= -1;
                for(int i=0;i<=n;i++){
                    row_aux[i] = A[ind_i][i];
                    A[ind_i][i] = A[q][i];
                    A[q][i] = row_aux[i];
                }
            }
            ////////////////////////////////////////////////

            // Swap cols
            double col_aux[n];
            if(q!=ind_j){
             nSwaps *= -1;
              for(int i=0;i<n;i++){
                    col_aux[i] = A[i][ind_j];
                    A[i][ind_j] = A[i][q];
                    A[i][q] = col_aux[i];
              }

             // Swap order of x in x_index as columns were swaped
             int x_aux = x_index[ind_j];
             x_index[ind_j] = x_index[q];
             x_index[q] = x_aux;
            }
         //////////////////////////////////////////////////


         // Convert to Upper triangular matrix
         for(int p=q+1;p<n;p++){
              for(int r=q+1;r<=n;r++){
                    A[p][r] = A[p][r]-(A[p-1-up][r]/A[q][q])*A[p][q];
                }
                A[p][q]=0;
                up++;
             }
        }
        //////////////////////////////////////////////////////////

        // Calculates the determinant
        det_A = nSwaps*determinant(A,n); //pow(-1.0,nSwaps)*

        // If the determinant could be calculated and it is different from 0 it means x could be calculated
        if(isnan(det_A) || det_A == 0){
            printf("Error solving x: check if matrix A is invertible or try rearranging its rows\n\n");

        }
        else
        {
            x = solver_TSupAug(A,n);
            // Rearrange x positions
            for(int i=0;i<n;i++){
                int index = x_index[i];
                x_sort[index] = x[i];
            }
            ////////////////////////////

            printf("\nThe values of the unknowns are:\n");
            for(int i=0;i<n;i++) printf("%lf\n",x_sort[i]);
            printf("\nThe determinant of the matrix A is:\n%lf\n\n",det_A);
        }
    }
    free(x);
    free(x_sort);
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
*       - long double* x: x solutions to augmented upper atrix A
*/
    double* x = malloc(sizeof(long double) * n);
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
