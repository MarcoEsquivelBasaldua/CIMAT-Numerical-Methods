Code name: Power_Iteration.c
Author: Marco Antonio Esquivel Basaldua
Date: 05/09/19

Code description:
- This Code captures Matrix A,
- captures the values of solutions to n linear equations in the matrix b,
- This Code captures Matrix A and finds its maximum eigenvalue and its eigenvector assotiated using the power iteration method
    
Inputs:
- M_BIG.txt
    File contenant the matrix from we will get the greatest eigenvalue and its eigenvector
    the first row in M_BIG.txt has two numbers referring to number of rows and columns in the matrix
    the matrix is then presented by rows and columns.

    
Outputs:
- None

Compiler instructions:
- By positioning in the folder where Power_Iteration.c and M_BIG.txt are stored, use gcc compiler: 
    gcc Power_Iteration.c -lm
- Run: 
    ./a.out M_BIG.txt
- Results:
    maximum eigenvalue of matrix entered
    eigenvector assotiated to eigenvalue found
    Results after execution are displayed in the terminal window.
