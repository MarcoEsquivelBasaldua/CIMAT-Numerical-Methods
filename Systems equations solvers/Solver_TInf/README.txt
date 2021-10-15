Code name: Solver_TInf.c
Author: Marco Antonio Esquivel Basaldua
Date: 25/08/19

Code description:
- This Code captures Matrix A,
- captures the values of solutions to n linear equations in the matrix b,
- If A is a square lower triangular Matrix and b has sufficient values as solutions for the linear system:
     get the determinant of the matrix A 
     calculates the values of x from the equation Ax = b
    
Inputs:
- M_TINF.txt
    File contenant the coefficients matrix from the equation Ax = b
    the first row in M_TINF.txt has two numbers referring to number of rows and columns in the coefficients matrix
    the coeficient matrix is then presented by rows and columns.
- V_TINF.txt
    File contenant the solutions b to the linear equations in Ax = b
    the first row in V_TINF.txt has two numbers referring to number of rows and columns in the matrix
    the matrix is then presented by rows and columns.
    
Outputs:
- None

Compiler instructions:
- By positioning in the folder where Solver_TInf.c, M_TINF.txt and V_TINF.txt are stored, use gcc compiler: 
    gcc Solver_TInf.c
- Run: 
    ./a.out M_TINF.txt V_TINF.txt
- Results:
    Determinant of matrix A
    x values in the equation Ax = b
    Results after execution are displayed in the terminal window.
