/*
*   Code name: quadratic_regression.c
*   Author: Marco Antonio Esquivel Basaldua
*   Date: 23/11/19
*
*   Code description:
*       This code captures the 101 values of x and y in a cloud of points,
*       then it computes and plots the fitting quadratic function corresponding
*       to least squares method
*/
#include <stdio.h>
#include <stdlib.h>

void quadratic_reg(int,double*,double*);       // Least squares function

int main( int argc , char* argv[] ){
    // Read .txt file to take interpolation points
    FILE* fin = NULL;
    // File with interpolation points is read and displayed
	fin = fopen(argv[1], "r");
	if(!fin) printf("ERROR: file not open %s\n", argv[1]);

    // Total data points to be read
    int n = 101;

    // Arrays with values of x and y
    double *x;
    x = (double*) malloc(n*sizeof(double));
    double *y;
    y = (double*) malloc(n*sizeof(double));

    for(int i=0; i<n; i++){
        fscanf(fin, "%lf", &x[i]);
        fscanf(fin, "%lf", &y[i]);

        //printf("%lf %lf\n",x[i],y[i]);
    }

    // Least squares method
    quadratic_reg(n,x,y);


    free(x);
    free(y);
    return 0;
}

void quadratic_reg(int n, double *x, double *y){
/*
*   Inputs:
*       - int n: total of points gotten
*       - double *x: array of x coordinates
*       - double *y: array of y coordinates
*
*   Outputs:
*       - none
*/

    double sum_x,sum_x2,sum_x3,sum_x4,sum_y,sum_xy,sum_x2y;
    sum_x = sum_x2 = sum_x3 = sum_x4 = sum_y = sum_xy = sum_x2y = 0;

    for(int i=0; i<n; i++){
        double xi = x[i];
        double yi = y[i];

        sum_x += xi;
        sum_x2 += (xi * xi);
        sum_x3 += (xi * xi * xi);
        sum_x4 += (xi * xi * xi * xi);
        sum_y += yi;
        sum_xy += (xi * yi);
        sum_x2y += (xi * xi * yi);
    }

    // a and b are computed from the function g(x) = a + bx + cx^2
    double c = ((sum_x2-sum_x*sum_x/n)*(sum_x2y-sum_x2*sum_y/n)-(sum_x3-sum_x2*sum_x/n)*(sum_xy-sum_x*sum_y/n))/((sum_x2-sum_x*sum_x/n)*(sum_x4-sum_x2*sum_x2/n)-(sum_x3-sum_x2*sum_x/n)*(sum_x3-sum_x2*sum_x/n));
    double b = ((sum_xy-sum_x*sum_y/n)*(sum_x4-sum_x2*sum_x2/n)-(sum_x2y-sum_x2*sum_y/n)*(sum_x3-sum_x2*sum_x/n))/((sum_x2-sum_x*sum_x/n)*(sum_x4-sum_x2*sum_x2/n)-(sum_x3-sum_x2*sum_x/n)*(sum_x3-sum_x2*sum_x/n));
    double a = (sum_y-b*sum_x-c*sum_x2)/n;

    printf("The fitting line is: g(x) = %lf + %lfx + %lfx^2\n",a,b,c);

    // Error
    double Error = 0;
    for(int i=0; i<n; i++)
        Error += (y[i]-(a+b*x[i]+c*x[i]*x[i]))*(y[i]-(a+b*x[i]+c*x[i]*x[i]));

    printf("The total Error is: %g\n",Error);

    // Create .txt file with points to graph
    char *nombre = "Exit points.txt";     // Name to exit file

    FILE* file = NULL;
	file = fopen( nombre, "w" );
	if(!file) printf("ERROR: file not open %s\n", nombre);

    for(int i=0; i<n; i++) fprintf(file, "%lf %lf\n",x[i],a+b*x[i]+c*x[i]*x[i]);
    fclose(file);
}
