/*
*   Code name: Max_Index.c
*   Author: Marco Antonio Esquivel Basaldua
*   Date: 21/08/19
*
*   Code description:
*   This Code captures the elements a[i][i] of the square Matrix A,
*   then calculates the maximum absolute value in the matrix
*   and stores its location
*/
#include <stdio.h>
#include <stdlib.h>

int* max_index(double**,int,int);   // Function to find maximum absolute value in A and returns its location

int main( int argc , char* argv[] ){
    int n,m;        // Size of matrix A
    double **A;     // Matrix A
    int i,j;        // Integers used in all for cycles
    int *index;     // Vector where solutions to unknowns x are stored

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

    index = max_index(A,n,m);
    printf("\nLocation of max abs value in A (rows,cols):(%d,%d)\n\n",index[0],index[1]);

    return 0;
}

int* max_index(double** A, int n, int m){
/*
*   Inputs:
*       - double **A: Pointer to the matrix A
*       - int n: number of rows in A
*       - int m: number of columns in A
*
*   Output:
*       - int* index: this array stores the location of the maximum number un matrix A
*/
    int* index = malloc(sizeof(int) * 2);
    double max = 0, aux = 0;
    int i,j;

    for(i=0;i<n;i++){
        for(j=0;j<m;j++){
            if(A[i][j] < 0) aux = -A[i][j];
            else aux = A[i][j];

            if(aux > max){
                max = aux;
                index[0] = i;
                index[1] = j;
            }
        }
    }
    i = index[0];
    j = index[1];

    // Prints the maximum absolute value and returns its location in the matrix A
    printf("The maximum absolute value is:\n%lf\n",A[i][j]);
    return index;
}
