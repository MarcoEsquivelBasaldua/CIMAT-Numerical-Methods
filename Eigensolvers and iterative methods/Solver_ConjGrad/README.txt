Code name: Solver_ConjGrad.c
Author: Marco Antonio Esquivel Basaldua
Date: 21/09/19

Code description:
This Code captures Matrix A,
captures the values of solutions to the n linear equations in the matrix b,
If A is a square Matrix and b has sufficient values as solutions for the linear system:
calculates the values of x from the equation Ax = b, by the conjugate gradient method that is an iterative method
used to find the values of the unknowns
    
Inputs:
- M_BIG.txt
    File contenant the coefficients matrix from the equation Ax = b
    the first row in M_BIG.txt has two numbers referring to number of rows and columns in the coefficients matrix
    the coeficient matrix is then presented by rows and columns.
- V_BIG.txt
    File contenant the solutions b to the linear equations in Ax = b
    the first row in V_BIG.txt has two numbers referring to number of rows and columns in the matrix
    the matrix is then presented by rows and columns.
    
Outputs:
- None

Compiler instructions:
- By positioning in the folder where Solver_ConjGrad.c, M_BIGL.txt and V_BIG.txt are stored, use gcc compiler: 
    gcc Solver_ConjGrad.c -lm
- Run: 
    ./a.out M_BIG.txt V_BIG.txt
- Results:
    x values in the equation Ax = b
    Results after execution are displayed in the terminal window.
