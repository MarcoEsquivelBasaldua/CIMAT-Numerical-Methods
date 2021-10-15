/*
*   Code name: Solver_Jacobi.c
*   Author: Marco Antonio Esquivel Basaldua
*   Date: 05/09/19
*
*   Code description:
*   This Code captures Matrix A,
*   captures the values of solutions to the n linear equations in the matrix b,
*   If A is a square Matrix and b has sufficient values as solutions for the linear system:
*      calculates the values of x from the equation Ax = b, by Jacobi's method that is an iterative method
*      used to find the values of the unknowns
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void Solver_Jacobi(double**, int, int, double**,int, int, double, int);      // Function to find the values of the unknowns by Jacobi's method
int is_possible(double**,int,int,double**,int,int);                         // Function to check if x could be found given A and b

int main( int argc , char* argv[] ){
    int n,m;                // Size of matrix A
    int k,l;                // Size of matrix b (solutions to linear equations)
    double **A;             // Matrix A
    double **b;             // Values of solutions to the n linear equations
    double toler = 1e-9;    // A tolerance will determine if the code has found the solutions
    int max_iter = 10000;   // This value is the maximum of iterations the code will run. This avoids the code to execute forever
    int i,j;                // Integers used in all for cycles

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

    for(i=0; i<k; i++) for(j=0; j<l; j++) fscanf(fin, "%lf", &b[i][j]);

    if( k<=20 && l<=20){                                             // Matrix b is only displayed if it is smaller than 20 rows and cols
        for(i=0; i<k; i++){
		    for(j=0; j<l; j++){
			    printf("%lf ", b[i][j]);
		    }
		    printf("\n");
	    }
    }
    else printf("Matrix b is too large to be displayed\n");
    printf("\n");
	fclose( fin );

    // Funtion to solve linear system
    Solver_Jacobi(A, n, m, b, k, l, toler, max_iter);

    // Memory is freed
    for(int i=0;i<n;i++) free(A[i]);
    free(A);
    for(int i=0;i<n;i++) free(b[i]);
    free(b);
    return 0;
}

void Solver_Jacobi(double** A, int n, int m, double** b, int k, int l, double toler, int max_iter){
/*
*   Inputs:
*       - double **A: Pointer to the matrix A
*       - double **b: Pointer to the matrix b
*       - int n: number of rows in A
*       - int m: number of columns in A
*       - int k: number of rows in b
*       - int l: number of columns in b
*       - double toler: tolerance near to cero to determine if the code has found the solutions
*       - int max_iter: This value is the maximum of iterations the code will run. This avoids the code to execute forever
*
*   Outputs:
*       - None
*/
    // First evaluates if this process can be done
    int flag;
    flag = is_possible(A,n,m,b,k,l);

    if(flag == 1){
        double* xO = malloc(sizeof(double) * n);    // Solutions to unknowns in the previous iteration are stored here
        double* xN = malloc(sizeof(double) * n);    // Solutions to unknowns in the current iteration are stored here

        // xO initialization, all entries are set to 1
        for(int i=0; i<n; i++) xO[i] = 1.0;

        int count_iter = 1;     // This variable keeps track on iterations currently taken
        double aux = 1000;      // Aux will be used to find the diference between xO and xN. This will help determine if the solutions were found

        // Iterations begin
        // while count_iter has not reached the maximum iterations permitted and aux is still greather than the permitted tolerance the cose keeps calculating
        while(count_iter <= max_iter && aux > toler){
            aux = 0;                                    // aux has to start at 0 at every iteration

            // New values for the unknowns are calculated given the values gotten in the previous iteration  according to Jaconi's method
            for(int i=0; i<n; i++){
                double sum =0;

                for(int j=0; j<n; j++) if(i != j) sum += A[i][j]*xO[j];

                xN[i] = (b[i][0]-sum)/A[i][i];

                if(xN[i] != 0) aux += (xN[i]-xO[i])*(xN[i]-xO[i])/(xN[i]*xN[i]);        // if condition in this line avoids division by cero
            }

            aux = sqrt(aux);

            // Evaluation to set xN as solution or keep iterations going
            if(aux > toler){
                if(count_iter < max_iter) for(int i=0; i<n; i++) xO[i] = xN[i];                     // xO is prepeared for next iteration taking xn values
                else printf("\nSolution does not converges\nIterations taken:%d\n\n",count_iter);   // iterations do not find the solution
            }
            else{
                printf("The values of the unknowns are:\n");                                        // xN is set as solution found
                for(int i=0; i<n; i++) printf("%.10lf\n",xN[i]);
                printf("\nIterations taken:%d\n",count_iter);
                printf("\nError:%.10lf\n\n",aux);
            }

            count_iter++;       // one more iteration is registered
        }

        free(xO);   // Memory is freed
        free(xN);
    }
}

int is_possible(double** A, int n, int m, double** b, int k, int l){
/*
*   Inputs:
*       - double **A: Pointer to the matrix A
*       - double **b: Pointer to the matrix b
*       - int n: number of rows in A
*       - int n: number of columns in A
*       - int n: number of rows in b
*       - int n: number of columns in b
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
