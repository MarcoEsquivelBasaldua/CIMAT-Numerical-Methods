/*
*   Code name: Jacobi.c
*   Author: Marco Antonio Esquivel Basaldua
*   Date: 14/09/19
*
*   Code description:
*   This Code captures Matrix A and finds its eigenvalues and their
*   eigenvectors assotiated using the Jacobi's method.
*
*   The results are stored in two .txt files:
*       "Eigenvalues.txt"
*       "Eigenvectors.txt"
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

void Jacobi(double**,double**,int,double,int);      // Jacobi method to find eigenvalues and eigenvectors
void max_index(double**,int,int,int*);             // Find biggest absolute value in a given matrix

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
    //int k = 10;                  // This is the number of eigenvalues we are interested in

    //double* lamda = malloc(sizeof(double) * k);        // All eigenvalues calculated will be here stored
    double **FI;                                       // Eigenvectors calculated will be here stored
    FI = calloc(n, sizeof *FI);
    for (int i=0; i<n; i++) FI[i] = calloc(n, sizeof *FI[i]);

    clock_t t_ini, t_fin;
    double secs;
    t_ini = clock();

    Jacobi(A,FI,n,toler,maxIter);

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

    fprintf(file,"Eigenvalues are:\n");
    for(int i=0; i<n; i++) fprintf(file,"%lf\n", A[i][i]);
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

    fprintf(file2,"Eigenvectors assotiated by cols:\n");
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++) fprintf(file2,"%lf   ",FI[i][j]);
        fprintf(file2,"\n");
    }
    fclose(file2);
    ////


    // Memory is freed
    for(int i=0;i<n;i++) free(A[i]);
    free(A);
    for(int i=0;i<n;i++) free(FI[i]);
    free(FI);

    return 0;
}

void max_index(double** A, int n, int m, int* index){
/*
*   Inputs:
*       - double **A: Pointer to the matrix A
*       - int n: number of rows in A
*       - int m: number of columns in A
*
*   Output:
*       - int* index: this array stores the location of the maximum number in matrix A, off the principal diagonal of A
*/
    double max = 0, aux = 0;
    int i,j;

    for(i=0;i<n;i++){
        for(j=0;j<m;j++){
            if(i != j){
                aux = fabs(A[i][j]);
                if(aux > max){
                    max = aux;
                    index[0] = i;
                    index[1] = j;
                }
            }
        }
    }
}

void Jacobi(double** A,double** FI,int n,double toler,int maxIter){
/*
*   Inputs:
*       - double **A: square matrix to calculate it eigenvalues and assotiated eigenvectors
*       - double **FI: mxn matrix where eigenvectors will be stored
*       - int n: dimension of square matrix A
*       - double toler: tolerance near to cero to determine if the code has found the solutions
*       - int max_iter: This value is the maximum of iterations the code will run. This avoids the code to execute forever
*
*   Outputs:
*       - double **A: eigenvalues calculated will be return by deferring
*       - double **FI: eigenvectors calculated will be return by deferring
*/

    // FI is initialized as the identity matrix
    for(int i=0;i<n;i++) FI[i][i] = 1;

    int *index;
    index = malloc(sizeof(int) * 2);

    for(int iter=1; iter<maxIter; iter++){
        //Find the position of the biggest value off the diagonal of A and store its location in the array "index"
        max_index(A,n,n,index);

        // i and j are the location of the biggest absolute value in A off diagonal
        int i = index[0];
        int j = index[1];

        // Test to see if the solution has converged
        if(fabs(A[i][j]) <= toler) break;
        else{
            double theta = 0.5*atan2( 2*A[i][j] , A[i][i]-A[j][j] );

            double sTheta = sin(theta);
            double cTheta = cos(theta);

            for(int l=0;l<n; l++){
                double tempA_i = A[l][i];
                double tempA_j = A[l][j];

                A[l][i] = tempA_i*cTheta + tempA_j*sTheta;
                A[l][j] = -tempA_i*sTheta + tempA_j*cTheta;

                double tempFI_i = FI[l][i];
                double tempFI_j = FI[l][j];

                FI[l][i] = tempFI_i*cTheta + tempFI_j*sTheta;
                FI[l][j] = -tempFI_i*sTheta + tempFI_j*cTheta;
            }

            for(int m=0;m<n;m++){
                double tempA_i = A[i][m];
                double tempA_j = A[j][m];

                A[i][m] = tempA_i*cTheta + tempA_j*sTheta;
                A[j][m] = -tempA_i*sTheta + tempA_j*cTheta;
            }
        }
    }

    free(index);
}
