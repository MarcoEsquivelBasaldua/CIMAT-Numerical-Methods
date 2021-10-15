Code name: Inverse_power.c
Author: Marco Antonio Esquivel Basaldua
Date: 14/09/19

Code description:
This Code captures Matrix A and finds its minimum eigenvalue and its eigenvector
assotiated using the Inverse power iteration method
    
Inputs:
- M_BIG.txt
    File contenant the matrix from we will get the smallest eigenvalue and its eigenvector
    the first row in M_BIG.txt has two numbers referring to number of rows and columns in the matrix
    the matrix is then presented by rows and columns.

    
Outputs:
- None

Compiler instructions:
- By positioning in the folder where Inverse_power.c and M_BIG.txt are stored, use gcc compiler: 
    gcc Inverse_power.c -lm
- Run: 
    ./a.out M_BIG.txt
- Results:
    minimum eigenvalue of matrix entered
    eigenvector assotiated to eigenvalue found
    Results after execution are displayed in the terminal window.
