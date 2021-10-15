/*
*   Code name: vectSum.c
*   Author: Marco Antonio Esquivel Basaldua
*   Date: 05/10/19
*
*   Code description:
*   This Code realizes the sum of two vectors and stores the result in another vector
*	using parallelization.
*	The procedure is done using 1, 2, 3 and 4 threads
*
*   The time results for each number of thread are stored in one .txt files ("Results.txt")
*   and ploted in one graph ("Sum of Vectors.png") for comparison.
*/

#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>
#include <math.h>

#define SIZE 25000000		// Number of elements on each vector

void suma(double* a, double* b, double*c, int size, int n_per_thread)
{
	int i;
	#pragma omp parallel for shared(a, b, c) private(i) schedule(dynamic, n_per_thread)
	for (i = 0; i < size; ++i)
	{
		c[i] = a[i] + b[i];
	}
}

double normV(double* V,int n){
/*
*   Intputs:
*       - double* V: vector to calculate it norm
*       - int n: number of entries in vector V
*
*   Outputs:
*       - double : euclidean norm of vector V
*/
    double norm = 0;
    for(int i=0; i<n; i++) norm += V[i]*V[i];

    norm = sqrt(norm);

    return norm;
}


int main(void){

	double* a = (double*)malloc(sizeof(double)*SIZE);
	double* b = (double*)malloc(sizeof(double)*SIZE);
	double* c = (double*)malloc(sizeof(double)*SIZE);

	clock_t t_ini, t_fin;
    double secs[4];					// Time taken for every number of thread is stored here
	int n_per_thread;				// determine how many elements each process will work on

	// Procedure is done with 1 to 4 threads
	for(int i=0; i<4; i++){
		omp_set_num_threads(i+1);
    	t_ini = clock();
		n_per_thread = SIZE/(i+1);
		suma(a, b, c, SIZE,n_per_thread);
		t_fin = clock();
    	secs[i] = (double)(t_fin - t_ini) / CLOCKS_PER_SEC;

		double norm = normV(c,SIZE);

	    printf("%.6lf seconds taken with %d threads\n", secs[i],i+1);
		printf("Norm is: %lf\n\n",norm);
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

	free(c);
	free(b);
	free(a);

	return 0;
}
