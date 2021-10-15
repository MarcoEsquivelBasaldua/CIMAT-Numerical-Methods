Code name: Solver_Doolittle.c
Author: Marco Antonio Esquivel Basaldua
Date: 25/08/19

Code description:
- This Code captures Matrix A,
- captures the values of solutions to n linear equations in the matrix b,
- If A is a square Matrix and b has sufficient values as solutions for the linear system:
     get the determinant of the matrix A 
     calculates the values of x from the equation Ax = b 
        - by first finding A = LU factorization using Doolittle algorithm
        - then finding solutions to y in the form Ly = b
        - finally finding solutions to x in the form Ux = y
    
Inputs:
- M_SMALL.txt
    File contenant the coefficients matrix from the equation Ax = b
    the first row in M_SMALL.txt has two numbers referring to number of rows and columns in the coefficients matrix
    the coeficient matrix is then presented by rows and columns.
- V_SMALL.txt
    File contenant the solutions b to the linear equations in Ax = b
    the first row in V_SMALL.txt has two numbers referring to number of rows and columns in the matrix
    the matrix is then presented by rows and columns.
    
Outputs:
- None

Compiler instructions:
- By positioning in the folder where Solver_Doolittle.c, M_SMALL.txt and V_SMALL.txt are stored, use gcc compiler: 
    gcc Solver_Doolittle.c
- Run: 
    ./a.out M_SMALL.txt V_SMALL.txt
- Results:
    Determinant of matrix A
    x values in the equation Ax = b
    Results after execution are displayed in the terminal window.
