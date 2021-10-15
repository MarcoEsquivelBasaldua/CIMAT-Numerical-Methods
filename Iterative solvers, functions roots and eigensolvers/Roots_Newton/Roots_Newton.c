/*
*   Author: Marco Antonio Esquivel Basaldua
*   Date: 28/08/19
*
*   Code description:
*       given a math function this code gives the location of every root in the funtcion
*       for this, users have to enter the number of roots to be calculated and give an interval on wich every root is located
*       this function uses Newton methos to find the roots
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double f_x(int,double);               // The mathematical function to work with is chosen and the value of x is given so its y coordinate is calculated
double f_x_prime(int,double);        // This function calculates the derivative of the math function selected inthe menu
void find_Roots(int,int,double);    // Function to locate roots by Newton-Raphson method

int main(void){
    int max_iter = 100;                      // Maximum number of trys on finding a root
    double toler = 0.0000000000001;         // Limit on the bigest value for the function to be said that a rooot was found

    int f_Selc;                                   // Function selected in the menu
    printf("FUNCTION MENU\n");
    printf("1) x^2          in the range [-10, -1]\n");
    printf("2) x^2-2        in the range [-10, 10]\n");
    printf("3) sin(x)       in the range [-pi, pi]\n");
    printf("4) 1/x^2        in the range [-1,1]\n");
    printf("5) x^3+3x^2+2x  in the range [-3,3]\n");
    printf("6) 1/x^2        in the range [-10000,10000]\n");

    printf("\nEnter the number of the function in the menu to work with:\n");

    scanf("%d",&f_Selc);

    // Funtion to find the roots
    find_Roots(f_Selc, max_iter, toler);

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
            return x*x;
            break;
        case 2:
            return x*x - 2;
            break;
        case 3:
            return sin(x);
            break;
        case 4:
            return 1/(x*x);
            break;
        case 5:
            return x*x*x + 3*x*x +2*x;
            break;
        case 6:
            return 1/(x*x);
            break;
        default:
            return NAN;
            break;
    }
}

double f_x_prime(int f_selec,double x){
/*
*   Inputs:
*       - int f_selec: this value corresponds to the math function to work with according to the menu
*       - double x: the value 'x' in the function displayed
*
*   Outputs:
*       - double : value for the function in the derivative for the given value 'x'
*/
    double delta = 0.0000001;
    return (f_x(f_selec,x+delta)-f_x(f_selec,x))/delta;
}

void find_Roots(int f_selec, int max_iter, double toler){
/*
*   Inputs:
*       - int f_selec: this value selcts tha math function to work with in the menu
*       - double* guess: each guess near the location of a root
*       - int max_iter: maximum number of trys on finding a root
*       - double toler: Limit on the bigest value for the function to be said that a rooot was found
*
*   Outputs:
*       None
*/
    if(f_selec >=1 && f_selec <=6){
        double xold = 0;              // A starting point approximation to start finding the root
        double root = NAN;                 // The found root will be stored here

        double l_min;              // Limits for the roots to be found are recupered from limits
        double l_max;

        switch (f_selec){         // minimum and maximum value in the window to start finding the rooot are initialized
        case 1:
            l_min = -10;
            l_max = -1;
            break;
        case 2:
            l_min = -10;
            l_max = 10;
            break;
        case 3:
            l_min = - 3.1415;
            l_max = 3.1415;
            break;
        case 4:
            l_min = -1;
            l_max = 1;
            break;
        case 5:
            l_min = -3;
            l_max = 3;
            break;
        case 6:
            l_min = -10000;
            l_max = 10000;
            break;
        default:
            break;
        }

        // Create file to graph the results

        // File with ponits to grph the function
        char *nombre = "f_x.txt";     // Name to exit file

        FILE* file = NULL;
	    file = fopen( nombre, "w" );
	    if(  !file  )  printf("Error: No se abrio %s\n" , nombre );

        int points = 160;
        double delta  = (l_max-l_min)/points;
        for(int i=0;i<points;i++) fprintf(file, "%lf %lf\n",l_min+i*delta, f_x(f_selec,l_min+i*delta));

        fclose(file);

        // Iterations begin to find the roots
        for(int iter = 1; iter <= max_iter; iter++){
            double xnew;

            xnew = xold - f_x(f_selec,xold)/f_x_prime(f_selec,xold);    // A new value for x is calculated by Newtons method

            double y_xnew = f_x(f_selec,xnew);                 // We get the y coordinate for the xnew value calculated


            if(fabs(y_xnew) <= toler){            // If this coordinate is equal to cero or suficiently close acordind to toler a root is found
                root = xnew;               // xnew is said to be the root found
                printf("\nRoot is located in:\n");
                printf("x = %lf\n", root);
                if(root < l_min || root > l_max) printf("But is not located in the display window.\n\n\n");
                break;
            }
            else if(iter == max_iter){          //If this method fails finding a root an error message is printed
                printf("\n Error: The root was not found\n");
            }
            else{
                xold = xnew;                   // Variables are prepared for next iteration
            }
        }

        // File with root points found
        nombre = "roots.txt";     // Name to exit file

        FILE* file2 = NULL;
	    file2 = fopen( nombre, "w" );
	    if(  !file2  )  printf("Error: No se abrio %s\n" , nombre );

        if(!isnan(root) && (root >= l_min && root <= l_max)) fprintf(file2, "%lf %lf",root, 0.0);

        fclose(file2);
    }
    else printf("The funtion you are trying is not in the menu\n\n"); // This warning appears if user does not type an option in the menu
}
