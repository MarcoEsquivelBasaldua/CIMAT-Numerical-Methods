/*
*   Code name: InterNewtonParticular.c
*   Author: Marco Antonio Esquivel Basaldua
*   Date: 17/10/19
*
*   Code description:
*       Given the function sin(x_i), this code gets four interpolation points as:
*       (x_0,sin(x_0)), (x_1,sin(x_1)), (x_1,sin(x_1)), (x_1,sin(x_1))
*       where:
*       x_0 = -3
*       x_ = -2.6
*       x_ = 0
*       x_ = 1.2
*
*       By using Newton method an interpolating polinomyal of degree less or equal to 3
*       is calculated and ploted in the same graph with
*       sin(x) and the interpolation points in the interval [-pi,pi]
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

void Newton(double*,int);       // Interpolation by Newton Method

int main(void){
    printf("\nFunction used to get interpolation points:\n");
    printf("3) sin(x)    in the range [-pi, pi]\n");

    int n = 4;                  // Total of interpolation points
    //double *x_i = malloc(n * sizeof(double));
    double x_i[] = {-3,-2.6,0,1.2};
    // Interpolation points
    printf("\nPoints used in interpolation:\n");
    for(int i=0; i<n; i++){
        printf("x%d = %lf   y%d = %lf\n",i,x_i[i],i,sin(x_i[i]));
    }

    printf("\nThe Interpolating Polynomial degree is less or equal to %d\n",n-1);

    // Interpolation by Newton Method
    // Time excecution starts
    clock_t t_ini, t_fin;
    double secs;
    t_ini = clock();

    Newton(x_i,n);
    printf("\n");

    // Time excecution is displayed
    t_fin = clock();
    secs = (double)(t_fin - t_ini) / CLOCKS_PER_SEC;
    printf("%.16g miliseconds taken in the algorithm\n\n", secs * 1000.0);

    return 0;
}

double nestedNewton(int i,double x,double **DD,double* x_i,int n){
    if(i == (n-1)) return DD[0][n-1];
    else return DD[0][i] + (x-x_i[i])*nestedNewton(i+1,x,DD,x_i,n);
}

void Newton(double* x_i,int n){
/*
*   Inputs:
*       - double* x_i: vector with interpolating points to be used
*       - int n: dimension of vector x_i
*
*   Outputs:
*       - none
*/
    double l_min = - 3.1415;              // Limits for the interval [a,b]
    double l_max = 3.1415;


    // Set of m points to be used in the graphs
    int m = 100;                             // number of points used to graph
    double x[m+1];                          // Values in the x axis
    double P[m+1];                         // Values in the y axis for P_n(x) (interpolated graph)
    double h = (l_max - l_min)/m;         // Distance between two consecutive x's

    // Values in the x axis to be used in the graphs
    for(int i=0; i<=m; i++){
        if(i == 0) x[i] = l_min;
        else x[i] = x[i-1] + h;
    }

    // Newton Method in action

    // Create matrix to store Divided Differences (DD)
    double **DD = (double **)malloc(n * sizeof(*DD));
    for(int j=0; j<n; j++) DD[j] = (double *) malloc((n-j) * sizeof(*DD[j]));

    for(int j=0; j<n; j++){
        int pos = j;
        for(int i=0; i<(n-j); i++){
            if(j == 0) DD[i][j] = sin(x_i[i]);
            else{
                DD[i][j] = (DD[i+1][j-1] - DD[i][j-1])/(x_i[pos] - x_i[i]);
            }
            pos++;
        }
    }

    /////
    printf("\nDivided Differences generated\n");
    for(int i=0; i<n; i++){
        for(int j=0; j<(n-i); j++){
            printf("%lf ",DD[i][j]);
        }
        printf("\n");
    }
///////

    for(int j=0; j<=m; j++){
        P[j] = nestedNewton(0,x[j],DD,x_i,n);
    }

    /* Error between f(x) and P(x)
    Calculated as the mean of the absolute values of the difference between
    every sample point in f(x) and P(x)*/
    double error = 0;
    for(int i=0; i<=m; i++){
        error += fabs(P[i] - sin(x[i]));
    }
    error /= m;
    printf("The error between the original function f(x) and the interpolating polinomyal P(x) is\n%lf\n",error);


    // Create .txt file with points to graph
    char *nombre = "Exit points.txt";     // Name to exit file

    FILE* file = NULL;
	file = fopen( nombre, "w" );
	if(!file) printf("ERROR: file not open %s\n", nombre);

    fprintf(file, "x, f(x), P_n(x)\n");
    for(int i=0; i<=m; i++) fprintf(file, "%lf %lf %lf\n",x[i], sin(x[i]),P[i]);
    fclose(file);

    // File with Interpolation points (commun points to both graphs)
    nombre = "Interpolation points.txt";     // Name to exit file

    FILE* file2 = NULL;
	file2 = fopen( nombre, "w" );
	if(!file2) printf("ERROR: file not open %s\n", nombre);

    fprintf(file, "x_i, f(x_i)\n");
    for(int i=0; i<n; i++) fprintf(file2, "%lf %lf\n",x_i[i], sin(x_i[i]));

    fclose(file2);

    for(int i=0;i<n;i++) free(DD[i]); // Se libera la memoria solicitada
    free(DD);
}
