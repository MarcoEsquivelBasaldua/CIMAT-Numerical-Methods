/*
*   Code name: Trapezoidal.c
*   Author: Marco Antonio Esquivel Basaldua
*   Date: 03/11/19
*
*   Code description:
*       User selects a function present in the start menu to get its integral value
*       by Simpson 3/8 method for the given interval.
*       Results are presented in the terminal window.
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

double f_x(int, double);                    // The mathematical function to work with is chosen and the value of x is given so its y coordinate is calculated
void Integral(int,int);                    // Function to calculate integral
double Simpson3_8(double,double,int,int); // Method used to get the integral


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

    // Specify the number of points to integrate the function
    printf("Enter the value of m:\n");
    int m;
    scanf("%d", &m);

    // n needs a multiple of 3, if it is not it must be adjusted
    if(m%3 != 0){
        if(m%3 == 1) m += 2;
        else if(m%3 == 2) m += 1;
        printf("n is adjusted to %d\n",m);
    }



    // Time excecution starts
    clock_t t_ini, t_fin;
    double secs;
    t_ini = clock();

    Integral(m,f_selc);
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

double Simpson3_8(double x_0,double h,int m,int f_selec){
/*
*   Inputs:
*       - double x_0: star of the interval to get the integral
*       - double h: distance between x_i and x_(i+1)
*       - int m: total of points to approximate the integral
*       - int f_selec: this value corresponds to the math function to work with according to the menu
*
*   Output:
*       - double integ: approximation to the integral by Simpson 3/8 method
*/
    double integ = 0;

    for(int i=0; i<=m; i++){
        double x_i = x_0 +i*h;
        if(i == 0 || i == m) integ += f_x(f_selec,x_i);
        else if(i%3 == 0) integ += 2*f_x(f_selec,x_i);
        else integ += 3*f_x(f_selec,x_i);
    }

    integ *= (3*h/8);

    return integ;
}

void Integral(int m, int f_selec){
/*
*   Inputs:
*       - int m: total of points to approximate the integral
*       - int f_selec: this value corresponds to the math function to work with according to the menu
*
*   Outputs:
*       - none
*/
    double x_0;        // Limits for the interval [x_0,x_m]
    double x_m;

    switch (f_selec){  // Minimum and maximum value in the interval to integrate according to function selected
    case 1:
        x_0 = -10;
        x_m = 10;
        break;
    case 2:
        x_0 = -10;
        x_m = 10;
        break;
    case 3:
        x_0 = - 3.14159;
        x_m = 3.14159;
        break;
    case 4:
        x_0 = -2;
        x_m = 4;
        break;
    case 5:
        x_0 = 0;
        x_m = 4*3.14159;
        break;
    case 6:
        x_0 = 1;
        x_m = 2;
        break;
    case 7:
        x_0 = 1;
        x_m = 5;
        break;
    case 8:
        x_0 = -1;
        x_m = 1;
        break;
    case 9:
        x_0 = 1;
        x_m = 2;
        break;
    default:
        break;
    }

    // Distance between x_i and x_(i+1)
    double h;
    if(m > 0) h =(x_m - x_0)/m;
    else h = 0;

    double integ;

    // Integral by Simson 3/8 rule
    integ = Simpson3_8(x_0,h,m,f_selec);

    printf("The integral is: %0.7g\n",integ);
}
