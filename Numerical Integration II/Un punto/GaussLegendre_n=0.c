/*
*   Code name: GaussLegendre_n=0.c
*   Author: Marco Antonio Esquivel Basaldua
*   Date: 05/11/19
*
*   Code description:
*       User selects a function present in the start menu to get its integral value
*       by Gauss-Legendre method for the given interval using 1 point (n=0).
*       Results are presented in the terminal window.
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define N 1

double f_x(int, double);                     // The mathematical function to work with is chosen and the value of x is given so its y coordinate is calculated
void Integral_n1(int);                     // Function to calculate integral
double Trapezoidal(double,double,int,int); // Method used to get the integral


int main(void){

    int f_selc;                                   // Function selected in the menu
    printf("\nFUNCTION MENU\n");
    printf("1) f(x) = 2 + x                in the range [-10, 10]\n");
    printf("2) f(x) = x^2-5                in the range [-10, 10]\n");
    printf("3) f(x) = sin(x)               in the range [-pi, pi]\n");
    printf("4) f(x) = x^2(5x-3)-2x^4+4x-5  in the range [-2,4]\n");
    printf("5) f(x) = sin(x + 1/2)         in the range [0,4pi]\n");
    printf("6) f(x) = e^(-x^2)             in the range [1,2]\n");
    printf("7) f(x) = x^3 + 2x^2           in the range [1,5]\n");
    printf("8) f(x) = 1/(1+x^2)            in the range [-1,1]\n");
    printf("9) f(x) = sqrt(1+cos²x)        in the range [1,2]\n");

    printf("\nEnter the number of the function in the menu to work with:\n");

    scanf("%d",&f_selc);
    while(f_selc < 1 && f_selc > 9){
        printf("ERROR: option %d is no avaiable, try again.\n",f_selc);
        scanf("%d",&f_selc);
    }


    // Time excecution starts
    clock_t t_ini, t_fin;
    double secs;
    t_ini = clock();

    Integral_n1(f_selc);
    printf("\n");

    // Time excecution is displayed
    t_fin = clock();
    secs = (double)(t_fin - t_ini) / CLOCKS_PER_SEC;
    printf("%lf miliseconds taken in the algorithm\n\n", secs * 1000.0);

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
        case 5:
            return sin(x + 0.5);
            break;
        case 6:
            return exp(-(x*x));
            break;
        case 7:
            return x*x*x + 2*x*x;
            break;
        case 8:
            return 1/(1+x*x);
            break;
        case 9:
            return sqrt(1+(cos(x) * cos(x)));
            break;
        default:
            return NAN;
            break;
    }
}

void Integral_n1(int f_selec){
/*
*   Inputs:
*       - int f_selec: this value corresponds to the math function to work with according to the menu
*
*   Outputs:
*       - none
*/
    double a;        // Limits for the interval [a,x_m]
    double b;

    switch (f_selec){  // Minimum and maximum value in the interval to integrate according to function selected
    case 1:
        a = -10;
        b = 10;
        break;
    case 2:
        a = -10;
        b = 10;
        break;
    case 3:
        a = - 3.14159;
        b = 3.14159;
        break;
    case 4:
        a = -2;
        b = 4;
        break;
    case 5:
        a = 0;
        b = 4*3.14159;
        break;
    case 6:
        a = 1;
        b = 2;
        break;
    case 7:
        a = 1;
        b = 5;
        break;
    case 8:
        a = -1;
        b = 1;
        break;
    case 9:
        a = 1;
        b = 2;
        break;
    default:
        break;
    }

    // x_i points
    double x_i[1] = {0};

    // Weights w_i
    double w_i[1] = {2};

    double integ = 0;

    for(int i=0; i<N; i++){
        integ += w_i[i] * f_x(f_selec,x_i[i]*((b-a)/2)+(a+b)/2);
    }

    integ *= ((b-a)/2);

    printf("The integral is: %0.7g\n",integ);
}
