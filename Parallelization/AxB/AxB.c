/*
*   Code name: AxB.c
*   Author: Marco Antonio Esquivel Basaldua
*   Date: 05/10/19
*
*   Code description:
*   This Code realizes the matrix matrix multiplication
*   and stores the result in one vector using parallelization.
*	The procedure is done using 1, 2, 3 and 4 threads
*
*   The time results for each number of thread are stored in one .txt files ("Results.txt")
*   and ploted in one graph ("Matrix times Matrix.png") for comparison.
*/

#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>
#include <string.h>

#define SIZE 800		// Number of rows and columns on each matrix

void AxB(double** A, double** B, double** C,int size, int n_per_thread){
    int i;
    #pragma omp parallel for shared(A, B, C) private(i) schedule(dynamic, n_per_thread)
    for(i=0; i<size; i++){
        for(int j=0; j<size; j++){
            double sum = 0;
            for(int k=0; k<size; k++) sum += A[i][k] * B[k][j];

            C[i][j] = sum;
        }
    }
}


int main(void){

	double ** A;
    A = (double**)malloc(SIZE * sizeof *A);
    for (int i=0; i<SIZE; i++) A[i] = (double*)malloc(SIZE * sizeof *A[i]);
    double ** B;
    B = (double**)malloc(SIZE * sizeof *B);
    for (int i=0; i<SIZE; i++) B[i] = (double*)malloc(SIZE * sizeof *B[i]);
    double ** C;
    C = (double**)malloc(SIZE * sizeof *C);
    for (int i=0; i<SIZE; i++) C[i] = (double*)malloc(SIZE * sizeof *C[i]);


	clock_t t_ini, t_fin;
    double secs[4];					// Time taken for every number of thread is stored here
	int n_per_thread;				// determine how many elements each process will work on

	// Procedure is done with 1 to 4 threads
	for(int i=0; i<4; i++){
		omp_set_num_threads(i+1);
    	t_ini = clock();
		n_per_thread = SIZE/(i+1);
		AxB(A, B, C, SIZE, n_per_thread);
		t_fin = clock();
    	secs[i] = (double)(t_fin - t_ini) / CLOCKS_PER_SEC;

	    printf("%.6lf seconds taken with %d threads\n", secs[i],i+1);
	}
	///////////////////////////////////////////////////////


	// Results are stored in .txt files
    printf("File with results is generated.\n\n");
    /////////
    char *nombre = "Results.txt";     // Name to out file

    FILE* file = NULL;
	file = fopen( nombre, "w" );
	if(  !file  ){
		printf("Error: No se abrio %s\n" , nombre );
	}

    fprintf(file,"Cores : Time\n");
    for(int i=0; i<4; i++)
		fprintf(file,"%d %g\n", i+1, secs[i]);
    fclose(file);
    ////

	for(int i=0;i<SIZE;i++) free(A[i]);
    free(A);

	return 0;
}
