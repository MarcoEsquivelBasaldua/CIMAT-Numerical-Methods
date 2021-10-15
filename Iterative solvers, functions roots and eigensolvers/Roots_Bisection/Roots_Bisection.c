/*
*   Code name: Roots_Bisection.c
*   Author: Marco Antonio Esquivel Basaldua
*   Date: 28/08/19
*
*   Code description:
*       given a math function this code gives the location of a root in the funtcion if it exists
*       for this, users have to decide which function to work with given a menu displayed in the terminal window
*       this method uses bisection to find the roots
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double f_x(int, double);             // The mathematical function to work with is chosen and the value of x is given so its y coordinate is calculated
void find_Roots(int,int,double);    // Function to locate roots by bisection

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

void find_Roots(int f_selec, int max_iter, double toler){
/*
*   Inputs:
*       - int f_selec: this value selcts tha math function to work with in the menu
*       - int max_iter: maximum number of trys on finding a root
*       - double toler: Limit on the bigest value for the function to be said that a rooot was found
*
*   Outputs:
*       None
*/
    double root = NAN;

    if(f_selec >=1 && f_selec <=6){
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
            l_min = - 3.1415;//M_1_PI;
            l_max = 3.1415;//M_1_PI;
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

        int points = 160;                           // This variable determines how many points will be graphed
        double delta  = (l_max-l_min)/points;
        for(int i=0;i<points;i++) fprintf(file, "%lf %lf\n",l_min+i*delta, f_x(f_selec,l_min+i*delta));

        fclose(file);

        // Iterations begin to find the roots
        for(int iter = 1; iter <= max_iter; iter++){
            double l_med = (l_max + l_min)/2;      // A new value between the limits is calculated

            double f_min = f_x(f_selec,l_min);             // y axis value is calculated for the given limits
            double f_med = f_x(f_selec,l_med);
            //double f_max = f_x(l_max);


            if(fabs(f_med) <= toler){       // If this coordinate is equal to cero or suficiently close acordind to toler a root is found
                root = l_med;         // l_med is said to be the root found
                printf("\nRoot is located in:\n");
                printf("x = %lf\n", root);
                break;
            }
            else if(iter == max_iter){      //If this method fails finding a root an error message is printed
                printf("\n Error: The root was not found\n");
            }
            else{
                if((f_min*f_med) > 0) l_min = l_med;        // Variables are prepared for next iteration
                else l_max = l_med;
            }

        }

        // File with root points found

        nombre = "roots.txt";     // Name to exit file

        FILE* file2 = NULL;
	    file2 = fopen( nombre, "w" );
	    if(  !file2  )  printf("Error: No se abrio %s\n" , nombre );

        if(!isnan(root)) fprintf(file2, "%lf %lf",root, 0);

        fclose(file2);

    }
    else printf("The funtion you are trying is not in the menu\n\n");  // This warning appears if user does not type an option in the menu

}
