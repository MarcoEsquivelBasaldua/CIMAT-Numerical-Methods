/*
*   Code name: Solver_QR.c
*   Author: Marco Antonio Esquivel Basaldua
*   Date: 21/09/19
*
*   Code description:
*   This Code captures Matrix A and finds its eigenvalues and eigenvectors
*   assotiated using the QR method iteration method
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

void Solver_QR(double**,double**,double**,double**,int,double,int);   // Function to find eigen values and eigenvectors
void QR(double**,double**,double**,int);                             // QR factorization
void max_index(double**,int,int,int*);                              // Find biggest absolute value in a given matrix off its diagonal
void matrix_multiply_mmd(double**,double**,double**,int,int,int);  // Function to multiply two matrices

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

    if( n<=20 && m<=20){                 // Matrix A is only displayed if it is smaller than 20 rows and cols
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

    double **Q;     // Matrix Q
    double **R;     // Matrix R

    Q = malloc(n * sizeof *Q);
    for (int i=0; i<n; i++) Q[i] = malloc(n * sizeof *Q[i]);

    R = calloc(n , sizeof *R);
    for (int i=0; i<n; i++) R[i] = calloc(n , sizeof *R[i]);

    double **FI;                                       // Eigenvectors calculated will be here stored
    FI = calloc(n, sizeof *FI);
    for (int i=0; i<n; i++) FI[i] = calloc(n, sizeof *FI[i]);

    double toler = 1e-9;                    // A tolerance will determine if the code has found the solutions
    int maxIter = 1000000;                   // This value is the maximum of iterations the code will run. This avoids the code to execute forever

    // Time excecution starts
    clock_t t_ini, t_fin;
    double secs;
    t_ini = clock();

    Solver_QR(A,Q,R,FI,n,toler,maxIter);
    printf("\n");

    // Time excecution is displayed
    t_fin = clock();
    secs = (double)(t_fin - t_ini) / CLOCKS_PER_SEC;
    printf("%.16g miliseconds taken in the algorithm\n", secs * 1000.0);

    // Results are stored in .txt files
    printf("File with Eigenvalues is generated.\n\n");
    /////////
    char *nombre = "Eigenvalues.txt";     // Name to out file

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
    char *nombre2 = "Eigenvectors.txt";     // Name to out file

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
    for(int i=0;i<n;i++) free(R[i]);
    free(R);
    for(int i=0;i<n;i++) free(Q[i]);
    free(Q);
    for(int i=0;i<n;i++) free(FI[i]);
    free(FI);

    return 0;
}

void Solver_QR(double** A, double** Q,double **R,double **FI, int n, double toler, int maxIter){
/*
*   Inputs:
*       - double **A: square matrix to calculate one of its eigenvalues and its assotiated eigenvector
*       - double **Q: Q factorization will be stored here
*       - double **R: R factorization will be stored here
*       - double **FI: FI eigenvectors will be stored here
*       - int n: number of rows in matrix A and vectors eVnew and eVold
*       - double toler: tolerance near to cero to determine if the code has found the solutions
*       - int max_iter: This value is the maximum of iterations the code will run. This avoids the code to execute forever
*
*   Outputs:
*       - double **A: eigenvalues will be stored in the diagonal of A
*       - double **FI: FI eigenvectors will be stored here returned by deferring
*/
    // FI is initialized as the identity matrix
    for(int i=0;i<n;i++) FI[i][i] = 1;

    double **FI_aux;                            // Matrix FI_aux
    FI_aux = malloc(n * sizeof *FI_aux);
    for (int i=0; i<n; i++) FI_aux[i] = malloc(n * sizeof *FI_aux[i]);

    int *index;
    index = malloc(sizeof(int) * 2);

    int iter;
    for(iter=1; iter<maxIter; iter++){
        QR(A,Q,R,n);
        matrix_multiply_mmd(R,Q,A,n,n,n);
        //Find the position of the biggest value off the diagonal of A and store its location in the array "index"
        max_index(A,n,n,index);

        matrix_multiply_mmd(FI,Q,FI_aux,n,n,n);

        // i and j are the location of the biggest absolute value in A off diagonal
        int i = index[0];
        int j = index[1];

        // Test to see if the solution has converged
        if(fabs(A[i][j]) <= toler){
            printf("Iterations taken: %d\n",iter);
            break;
        }
        else{
            for(int i=0;i<n;i++) for(int j=0; j<n; j++) FI[i][j] = FI_aux[i][j];

        }
    }

    free(index);
    for(int i=0;i<n;i++) free(FI_aux[i]);
    free(FI_aux);
}

void QR(double** A, double** Q,double **R, int n){
/*
*   Inputs:
*       - double** A: Matrix to be factored
*       - double **Q: Q factorization will be stored here
*       - double **R: R factorization will be stored here
*       - int n: dimension of A,Q and R
*/
    double sum;
    for(int k=0;k<n; k++){
        sum =0;
        for(int j=0; j<n; j++) sum += A[j][k]*A[j][k];
        R[k][k] = sqrt(sum);
        for(int j=0; j<n; j++) Q[j][k] = A[j][k]/R[k][k];
        for(int i=k+1; i<n; i++){
            sum = 0;
            for(int j=0; j<n; j++) sum += A[j][i]*Q[j][k];
            R[k][i] = sum;
            for(int j=0; j<n; j++) A[j][i] = A[j][i] - R[k][i]*Q[j][k];
        }
    }

    // Active to display Q and R
    /*for(int i=0;i<n;i++){
        for(int j=0; j<n; j++){
            printf("%lf ", R[i][j]);
        }
        printf("\n");
    }

    printf("\n");
    for(int i=0;i<n;i++){
        for(int j=0; j<n; j++){
            printf("%lf ", Q[i][j]);
        }
        printf("\n");
    }*/
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



/* Multiply matrix a x b = c and save the result in the matrix c */
void matrix_multiply_mmd(double ** a, double ** b, double ** c, int arows, int acols, int bcols){
    double ** a_i = a;
    double ** b_i = b;
    double ** c_i = c;

    for (int i = 0; i < arows; i++, a_i++, b_i++, c_i++) {
        double ** b_k = b;
        double * a_ik = *a_i;
        memset(*c_i, 0, acols * sizeof(double));
        for (int k = 0; k < acols; k++, a_ik++, b_k++) {
            double * b_kj = *b_k;
            double * c_ij = *c_i;
            for (int j = 0; j < bcols; j++, b_kj++, c_ij++) {
                *c_ij += (*a_ik) * (*b_kj);
            }
        }
    }
}
