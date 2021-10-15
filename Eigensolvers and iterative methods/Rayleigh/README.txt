Code name: Rayleigh.c
Author: Marco Antonio Esquivel Basaldua
Date: 14/09/19

Code description:
This Code captures Matrix A and finds one of its eigenvalues and its eigenvector
assotiated using the Rayleigh method

The results are displayed in the terminal screen
    
Inputs:
- M_BIG.txt
    File contenant the matrix from we will get the eigenvalue and its eigenvector
    the first row in M_BIG.txt has two numbers referring to number of rows and columns in the matrix
    the matrix is then presented by rows and columns.

    
Outputs:
- None

Compiler instructions:
- By positioning in the folder where Rayleigh.cc and M_BIG.txt are stored, use gcc compiler: 
    gcc Rayleigh.c -lm
- Run: 
    ./a.out M_BIG.txt
- Results:
    one of the eigenvalues of matrix entered
    eigenvector assotiated to eigenvalue found
