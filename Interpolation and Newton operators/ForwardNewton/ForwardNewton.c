/*
*   Code name: ForwardNewton.c
*   Author: Marco Antonio Esquivel Basaldua
*   Date: 26/10/19
*
*   Code description:
*       This code captures the value n that correspondes to the total amount of interpolation points (n+1),
*       next it captues the x_0 value and x_n values in the file"Enter points.txt".
*       According to user selection, points in "Enter points.txt" are interpolated in a given funtion by
*       Gregory-Newton forward method.
*       The original function, the interpolated function and the interpolation points are ploted on the same graph,
*       data to plot are stored in exit .txt files.
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

double f_x(int, double);                     // The mathematical function to work with is chosen and the value of x is given so its y coordinate is calculated
int factorial(int);                         // Factorial function
void ForwardNewton(double,double,int,int); // Interpolation by Gregory-Newton Forward method

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

    int n, x_0, x_n;                  // Start point, end point and number of points
    fscanf(fin, "%d", &n);
    fscanf(fin, "%d", &x_0);
    fscanf(fin, "%d", &x_n);

    fclose(fin);

    printf("\nThe Interpolating Polynomial degree is less or equal to %d\n\n",n-1);

    // Interpolation by Newton Method
    // Time excecution starts
    clock_t t_ini, t_fin;
    double secs;
    t_ini = clock();

    ForwardNewton(x_0,x_n,n,f_selc);
    printf("\n");

    // Time excecution is displayed
    t_fin = clock();
    secs = (double)(t_fin - t_ini) / CLOCKS_PER_SEC;
    printf("%.16g miliseconds taken in the algorithm\n\n", secs * 1000.0);

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

int factorial(int n){
/*
*   Inputs:
*       - int n: Integer number to calculate its factorial
*
*   Outputs:
*       - int: Factorial of n
*/
    int fact = 1;

    for (int i=1; i <= n; i++)
        fact = fact * i;

    return fact;
}

void ForwardNewton(double x_0, double x_n, int n, int f_selec){
/*
*   Inputs:
*       - double x_0: x value of first interpolation point
*       - double x_n: x value of last interpolation point
*       - int n: total of interpolation points
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

    // Distance between x_i and x_(i+1)
    double h;
    if(n != 0) h =(x_n - x_0)/(n-1);
    else h = 0;

    // Vector with interpolating points x_i
    double x_i[n];
    x_i[0] = x_0;
    for(int i=1; i<n; i++) x_i[i] = x_i[i-1] + h;

    // Set of m points to be used in the graphs
    int m = 100;                             // number of points used to graph
    double x[m+1];                          // Values in the x axis
    double P[m+1];                         // Values in the y axis for P_n(x) (interpolated graph)
    double h_plot = (b - a)/m;            // Distance between two consecutive x's

    // Values in the x axis to be used in the graphs
    for(int i=0; i<=m; i++){
        if(i == 0) x[i] = a;
        else x[i] = x[i-1] + h_plot;
    }

    // Gregory-Newton Forward Method in action

    // Create matrix to store forward differences
    double **FD = (double **)malloc(n * sizeof(*FD));
    for(int j=0; j<n; j++) FD[j] = (double *) malloc((n-j) * sizeof(*FD[j]));

    for(int j=0; j<n; j++){

        for(int i=0; i<(n-j); i++){
            if(j == 0) FD[i][j] = f_x(f_selec,x_i[i]);
            else FD[i][j] = FD[i+1][j-1] - FD[i][j-1];
        }
    }

    /////
    printf("\nForward Differences generated\n");
    for(int i=0; i<n; i++){
        for(int j=0; j<(n-i); j++){
            printf("%lf ",FD[i][j]);
        }
        printf("\n");
    }
///////

    double aux_2_old = 1;
    for(int j=0; j<=m; j++){
        double aux = 0;
        double aux_1;
        double aux_2;
        double aux_2_old = 1;
        for(int i=1; i<n; i++){
            aux_1 = FD[0][i] /(factorial(i)*pow(h,i));

            aux_2 = aux_2_old * (x[j] - x_i[i-1]);
            aux_2_old = aux_2;

            aux += aux_1 * aux_2;
        }
        P[j] = aux + FD[0][0];
    }

    /* Error between f(x) and P(x)
    Calculated as the mean of the absolute values of the difference between
    every sample point in f(x) and P(x)*/
    double error = 0;
    for(int i=0; i<=m; i++){
        error += fabs(P[i] - f_x(f_selec,x[i]));
    }
    error /= m;
    printf("\nThe error between the original function f(x) and the interpolating polinomyal P(x) is\n%lf\n",error);


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

    for(int i=0;i<n;i++) free(FD[i]); // Se libera la memoria solicitada
    free(FD);
}
