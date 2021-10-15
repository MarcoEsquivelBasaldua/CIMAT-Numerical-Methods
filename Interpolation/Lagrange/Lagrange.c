/*
*   Code name: Lagrange.c
*   Author: Marco Antonio Esquivel Basaldua
*   Date: 15/10/19
*
*   Code description:
*       This code captures n x position values in "Enter points.txt" corresponding to the interpolation points to work with.
*       According to user selection, points in "Enter points.txt" are interpolated in a given funtion by Lagrange method.
*       The original function, the interpolated function and the interpolation points are ploted on the same graph,
*       data to plot are stored in exit .txt files.
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

double f_x(int, double);             // The mathematical function to work with is chosen and the value of x is given so its y coordinate is calculated
void Lagrange(double*,int,int);     // Interpolation by Lagrange Method


int main(int argc, char* argv[]){

    int f_selc;                                   // Function selected in the menu
    printf("\nFUNCTION MENU\n");
    printf("1) f(x) = 2 + x                in the range [-10, 10]\n");
    printf("2) f(x) = x^2-5                in the range [-10, 10]\n");
    printf("3) f(x) = sin(x)               in the range [-pi, pi]\n");
    printf("4) f(x) = x^2(5x-3)-2x^4+4x-5  in the range [-2,4]\n");

    printf("\nEnter the number of the function in the menu to work with:\n");

    scanf("%d",&f_selc);
    while(f_selc < 1 && f_selc > 4){
        printf("ERROR: option %d is no avaiable, try again.\n",f_selc);
        scanf("%d",&f_selc);
    }


    // Read .txt file to take interpolation points
    FILE* fin = NULL;
    // File with interpolation points is read and displayed
	fin = fopen(argv[1], "r");
	if(!fin) printf("ERROR: file not open %s\n", argv[1]);

    int n;                  // Total of interpolation points
    fscanf(fin, "%d", &n);

    double *x_i = malloc(n * sizeof(double));
    // Points read
    printf("\nPoints used in interpolation:\n");
    for(int i=0; i<n; i++){
        fscanf(fin, "%lf", &x_i[i]);
        printf("x%d = %lf   y%d = %lf\n",i,x_i[i],i,f_x(f_selc,x_i[i]));
    }
    fclose(fin);

    printf("\nThe Interpolating Polynomial degree is less or equal to %d\n\n",n-1);

    // Interpolation by Lagrange Method

    // Time excecution starts
    clock_t t_ini, t_fin;
    double secs;
    t_ini = clock();

    Lagrange(x_i,n,f_selc);
    printf("\n");

    // Time excecution is displayed
    t_fin = clock();
    secs = (double)(t_fin - t_ini) / CLOCKS_PER_SEC;
    printf("%.16g miliseconds taken in the algorithm\n\n", secs * 1000.0);

    free(x_i);

    return 0;
}

double f_x(int f_selec, double x){
/*
*   Inputs:
*       - int f_selec: this value corresponds to the math function to work with according to the menu
*       - double x: the value 'x' in the function displayed
*
*   Outputs:
*       - double* :
*/
    switch (f_selec){
        case 1:
            return 2 + x;
            break;
        case 2:
            return x*x-5;
            break;
        case 3:
            return sin(x);
            break;
        case 4:
            return x*x*(5*x-3)-2*x*x*x*x+4*x-5;
            break;
        default:
            return NAN;
            break;
    }
}

void Lagrange(double* x_i,int n, int f_selec){
/*
*   Inputs:
*       - double* x_i: vector with interpolating points to be used
*       - int n: dimension of vector x_i
*       - int f_selec: this value corresponds to the math function to work with according to the menu
*
*   Outputs:
*       - none
*/

    double a;              // Limits for the interval [a,b]
    double b;

    switch (f_selec){         // minimum and maximum value in the window to plot according to user selection
    case 1:
        a = -10;
        b = 10;
        break;
    case 2:
        a = -10;
        b = 10;
        break;
    case 3:
        a = - 3.1415;
        b = 3.1415;
        break;
    case 4:
        a = -2;
        b = 4;
        break;
    default:
        break;
    }

    // Set of m points to be used in the graphs
    int m = 100;                             // number of points used to graph
    double x[m+1];                          // Values in the x axis
    double P[m+1];                         // Values in the y axis for P_n(x) (interpolated graph)
    double h = (b - a)/m;         // Distance between two consecutive x's

    // Values in the x axis to be used in the graphs
    for(int i=0; i<=m; i++){
        if(i == 0) x[i] = a;
        else x[i] = x[i-1] + h;
    }

    // Lagrange Method in action
    for(int j=0; j<=m; j++){

        double aux_y = 0;
        for(int i=0; i<n; i++){

            double li = 1;           // Lagrange coefficients
            for(int k=0; k<n; k++){
                if(k != i) li *= (x[j] - x_i[k])/(x_i[i] - x_i[k]);
            }

            aux_y += f_x(f_selec,x_i[i]) * li;
        }
        P[j] = aux_y;
    }

    /* Error between f(x) and P(x)
    Calculated as the mean of the absolute values of the difference between
    every sample point in f(x) and P(x)*/
    double error = 0;
    for(int i=0; i<=m; i++){
        error += fabs(P[i] - f_x(f_selec,x[i]));
    }
    error /= m;
    printf("The error between the original function f(x) and the interpolating polinomyal P(x) is\n%lf\n",error);


    // Create .txt file with points to graph
    char *nombre = "Exit points.txt";     // Name to exit file

    FILE* file = NULL;
	file = fopen( nombre, "w" );
	if(!file) printf("ERROR: file not open %s\n", nombre);

    fprintf(file, "x, f(x), P_n(x)\n");
    for(int i=0; i<=m; i++) fprintf(file, "%lf %lf %lf\n",x[i], f_x(f_selec,x[i]),P[i]);
    fclose(file);

    // File with Interpolation points (commun points to both graphs)
    nombre = "Interpolation points.txt";     // Name to exit file

    FILE* file2 = NULL;
	file2 = fopen( nombre, "w" );
	if(!file2) printf("ERROR: file not open %s\n", nombre);

    fprintf(file, "x_i, f(x_i)\n");
    for(int i=0; i<n; i++) fprintf(file2, "%lf %lf\n",x_i[i], f_x(f_selec,x_i[i]));

    fclose(file2);


}
