/*
*   Author: Marco Antonio Esquivel Basaldua
*   Date: 28/08/19
*
*   Code description:
*       given a math function this code gives the location of every root in the funtcion
*       for this, users have to enter the number of roots to be calculated and give an interval on wich every root is located
*       this method uses bisection to find the roots
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double f_x(int, double);                     // The mathematical function to work with is chosen and the value of x is given so its y coordinate is calculated
void find_Roots(int,double*,int,double);    // Function to locate roots by bisection

int main(void){
    int max_iter = 100;                      // Maximum number of trys on finding a root
    double toler = 0.0000000000001;         // Limit on the bigest value for the function to be said that a rooot was found

    int f_Selc;                                   // Function selected in the menu
    printf("Enter the function in te menu to work with:\n");

    scanf("%d",&f_Selc);

//////////////////////////
    double* limits = malloc(sizeof(double) * (n*2));        // limits stores an interval per root where each root is located
    printf("Enter a star limit and a end limit for each root to be located in ascendant order:\n");
    for(int i=0; i<n*2; i++) scanf("%lf", &limits[i]);
//////////////////////

    find_Roots(f_Selc, limits, max_iter, toler);

    free(limits);
    return 0;
}

double* f_x(int f_selec, double x){
/*
*   Inputs:
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

void find_Roots(int n, double* limits, int max_iter, double toler){
/*
*   Inputs:
*       - int n: number of roots to be located
*       - double* limits: each consecutive pair are the limits between every root is located
*       - int max_iter: maximum number of trys on finding a root
*       - double toler: Limit on the bigest value for the function to be said that a rooot was found
*
*   Outputs:
*       None
*/
    double* roots = malloc(sizeof(double) * n);     // Roots will be stored here

    for(int i=1; i <= n; i++){
        double l_min = limits[(i-1)*2];             // Limits for the roots to be found are recupered from limits
        double l_max = limits[(i-1)*2+1];

        for(int iter = 1; iter <= max_iter; iter++){
            double l_med = (l_max + l_min)/2;      // A new value between the limits is calculated

            double f_min = f_x(l_min);             // y axis value is calculated for the given limits
            double f_med = f_x(l_med);
            //double f_max = f_x(l_max);


            if(fabs(f_med) <= toler){       // If this coordinate is equal to cero or suficiently close acordind to toler a root is found
                roots[i-1] = l_med;         // l_med is said to be the root found
                break;
            }
            else if(iter == max_iter){      //If this method fails finding a root an error message is printed
                roots[i-1] = NAN;
                printf("\n Error: The root number %d was not found\n",i);
            }
            else{
                if((f_min*f_med) > 0) l_min = l_med;        // Variables are prepared for next iteration
                else l_max = l_med;
            }

        }
    }

    printf("\nRoots are located in:\n");                // Roots found are printed
    for(int i=0; i<n; i++) printf("x_%d = %lf\n", i+1, roots[i]);

    free(roots);
}
