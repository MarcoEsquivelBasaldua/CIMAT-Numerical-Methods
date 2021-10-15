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

double f_x(double);                          // This is the function whose roots will be located
void find_Roots(int,double*,int,double);    // Function to locate roots by Newton

int main(void){
    int max_iter = 100;                      // Maximum number of trys on finding a root
    double toler = 0.0000000000001;         // Limit on the bigest value for the function to be said that a rooot was found

    int n;                                   // Number of roots to be located
    printf("Enter the number of roots to be located:\n");
    scanf("%d",&n);

    double* guess = malloc(sizeof(double) * n);        // limits stores an interval per root where each root is located
    printf("Enter a good gues for each root to be found, the all guesses by increasing order:\n");
    for(int i=0; i<n; i++) scanf("%lf", &guess[i]);

    find_Roots(n, guess, max_iter, toler);

    free(guess);
    return 0;
}

double f_x(double x){
/*
*   Inputs:
*       - double x: the value 'x' in the function displayed
*
*   Outputs:
*       - double : value for the function in the 'y' axis for the given value 'x'
*/
    return cos(x);
}

double f_x_prime(double x){
/*
*   Inputs:
*       - double x: the value 'x' in the function displayed
*
*   Outputs:
*       - double : value for the function in the derivative for the given value 'x'
*/
    double delta = 0.0000001;
    return (f_x(x+delta)-f_x(x))/delta;
}

void find_Roots(int n, double* guess, int max_iter, double toler){
/*
*   Inputs:
*       - int n: number of roots to be located
*       - double* guess: each guess near the location of a root
*       - int max_iter: maximum number of trys on finding a root
*       - double toler: Limit on the bigest value for the function to be said that a rooot was found
*
*   Outputs:
*       None
*/
    double* roots = malloc(sizeof(double) * n);     // Roots will be stored here

    for(int i=1; i <= n; i++){
        double xold = guess[i-1];                      // For each root to be found we have to start by a guess given by the user

        for(int iter = 1; iter <= max_iter; iter++){
            double xnew;

            xnew = xold - f_x(xold)/f_x_prime(xold);    // A new value for x is calculated by Newtons method

            double y_xnew = f_x(xnew);                 // We get the y coordinate for the xnew value calculated


            if(fabs(y_xnew) <= toler){            // If this coordinate is equal to cero or suficiently close acordind to toler a root is found
                roots[i-1] = xnew;               // xnew is said to be the root found
                break;
            }
            else if(iter == max_iter){          //If this method fails finding a root an error message is printed
                roots[i-1] = NAN;
                printf("\n Error: The root number %d was not found\n",i);
            }
            else{
                xold = xnew;                   // Variables are prepared for next iteration
            }

        }
    }

    printf("\nRoots are located in:\n");        // Roots found are printed
    for(int i=0; i<n; i++) printf("x_%d = %lf\n", i+1, roots[i]);

    free(roots);
}
