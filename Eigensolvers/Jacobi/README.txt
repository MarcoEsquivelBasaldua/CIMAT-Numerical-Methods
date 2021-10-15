Code name: Jacobi.c
Author: Marco Antonio Esquivel Basaldua
Date: 14/09/19

Code description:
This Code captures Matrix A and finds its eigenvalues and their
eigenvectors assotiated using the Jacobi's method.

The results are stored in two .txt files:
       "Eigenvalues.txt"
       "Eigenvectors.txt"
    
Inputs:
- M_BIG.txt
    File contenant the matrix from we will get the eigenvalues and their eigenvectors
    the first row in M_BIG.txt has two numbers referring to number of rows and columns in the matrix
    the matrix is then presented by rows and columns.

    
Outputs:
- "Eigenvalues.txt"
- "Eigenvectors.txt"

Compiler instructions:
- By positioning in the folder where Jacobi.c and M_BIG.txt are stored, use gcc compiler: 
    gcc Jacobi.c -lm
- Run: 
    ./a.out M_BIG.txt
- Results:
    eigenvalues of matrix entered
    eigenvectors assotiated to eigenvalues found
    Results after execution are stored in .txt files.
