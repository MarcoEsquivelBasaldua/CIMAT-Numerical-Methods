Code name: Inverse_power_kSmallest.c
Author: Marco Antonio Esquivel Basaldua
Date: 14/09/19

Code description:
This Code captures Matrix A and finds its K smallest eigenvalues and their
eigenvectors assotiated using the inverse power iteration method.
The value K is defined in the code. To change this value go to line 19 in the code

The results are stored in two .txt files:
       "Eigenvalues.txt"
       "Eigenvectors.txt"
    
Inputs:
- M_BIG.txt
    File contenant the matrix from we will get the K smallest eigenvalues and their eigenvectors
    the first row in M_BIG.txt has two numbers referring to number of rows and columns in the matrix
    the matrix is then presented by rows and columns.

    
Outputs:
- "Eigenvalues.txt"
- "Eigenvectors.txt"

Compiler instructions:
- By positioning in the folder where Inverse_power_kSmallest.c and M_BIG.txt are stored, use gcc compiler: 
    gcc Inverse_power_kSmallest.c -lm
- Run: 
    ./a.out M_BIG.txt
- Results:
    K smallest eigenvalues of matrix entered
    eigenvectors assotiated to eigenvalues found
    Results after execution are stored in .txt files.
