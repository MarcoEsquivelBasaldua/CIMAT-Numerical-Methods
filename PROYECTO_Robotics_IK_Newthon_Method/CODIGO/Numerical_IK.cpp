/*
*   Code name: Numerical_IK.cpp
*   Author: Marco Antonio Esquivel Basaldua
*   Date: 18/12/19
*
*   Code description:
*       Given a RRR robotic arm in three dimensions which lenght of each link is specified in the
*       file "lengths.txt", this code determines the set of nedded angles to each link so the end
*       effector of the robot will pass throught the points specified in the file "Path.txt".
*       The angles generated can be used to simulate the robotic arm in the software MATLAB.
*/
#include <iostream>
#include <map>
#include <cmath>
#include <fstream>
#include <sstream>
#include <vector>
#include <ctime>

using namespace std;

typedef vector<vector<double> > Matrix;
#define PI 3.14159265358979323846

// New functions to this priject
void readLenghts(vector<double>&);                   // Read the dimensions from each link
void readPath(Matrix&);                             // Read the points the robot has to visit
void writeAngles(Matrix);                          // Write generated angles in .txt file
Matrix DH(double,double,double,double);           // Get 4x4 matrix using DH convention
Matrix Jacobian(vector<double>,vector<double>);  // Get 3x3 Jacobian matrix
Matrix IK(Matrix&,vector<double>&);             // Permorms Inverse Kinematics

// Already existant funtions
vector<double> matXvect(Matrix,vector<double>);               // Function to multiply a matrix times a vector
Matrix AXB(Matrix,Matrix);                                   // Function to multiply two matrices
vector<double> vectSubs(vector<double>, vector<double>);    // Function to substract two vectors
vector<double> vectAdd(vector<double>, vector<double>);    // Funtion to add two vectors
double normV(vector<double>);                             // Get euclidean norm from given vector
Matrix Inverse(Matrix,int, int);                         // Function to get inverse Matrix of A
void solver_TSupAug(Matrix&,int,vector<double>&);       // Function to solve x given an upper triangular augmented matrix
void solver_TInfAug(Matrix&,int,vector<double>&);      // Function to solve x given a lower triangular augmented matrix



int main(void){
    vector<double> lenghts;     // Stores the lenght of each link
    Matrix Path, angles;        // The points to be visited will be stored here

    readLenghts(lenghts);
    readPath(Path);

    angles = IK(Path,lenghts);
    writeAngles(angles);

    cout << "Angles generated and stored." << endl;

    return 0;
}

void readLenghts(vector<double> &lenghts){
/*
*   Inputs:
*       - vector<double> &lenghts: vector where lenght o each from three links will be stored
*
*   Outputs:
*       - vector<double> &lenghts: results are stored by rediffering
*/
    double l1,l2,l3;

    ifstream input;
    input.open("lengths.txt");
    if(input.is_open())
        input >> l1 >> l2 >> l3;

    lenghts.push_back(l1);
    lenghts.push_back(l2);
    lenghts.push_back(l3);
}

void readPath(Matrix& Path){
/*
*   Input:
*       - Matrix Path: matrix to store the points the robot will pass trought
*
*   Output:
*       - Matrix Path: results stored by rediffering
*/
    int m;
    Matrix vertices;
    vector<double> row;
    double val;

    row.push_back(0);
    row.push_back(-2);
    row.push_back(5);
    vertices.push_back(row);
    row.clear();

    ifstream input;
    input.open("Path.txt");
    if(input.is_open()){
        input >> m;
        for(int i=0; i<m; i++){
            for(int j=0; j<3; j++){
                input >> val;
                row.push_back(val);
            }
            vertices.push_back(row);
            row.clear();
        }
    }
    //vertices.push_back(vertices[0]);

    int n = 10;
    for(int i=0; i<m; i++){
        vector<double> u,v,w,p;
        u = vertices[i];
        v = vertices[i+1];
        w = vectSubs(u,v);
        double mag_w = sqrt(w[0]*w[0] + w[1]*w[1] + w[2]*w[2]);
        for(int k=0; k<3; k++)
            w[k] = w[k]/mag_w;

        double delta = mag_w/n;
        Path.push_back(u);
        for(int j=1; j<=n; j++){
            p.clear();
            for(int k=0; k<3; k++){
                p.push_back(u[k] - j*delta*w[k]);
            }
            Path.push_back(p);
        }
    }
}

void writeAngles(Matrix angles){
/*
*   Input:
*       - Matrix angles: set of nx3 angles for each link
*
*   Output:
*       - None
*/
    ofstream output;
    // Nombre del archivo
    output.open("angles.txt");

    for(int i=0; i<angles.size(); i++){
        for(int j=0; j<angles[i].size(); j++)
            output << angles[i][j] << " ";
        output << endl;
    }

    output.close();
}

Matrix DH(double a,double alpha,double d,double theta){
/*
*   Inputs
*       - double a,double alpha,double d,double theta
*
*   Outputs:
*       - Matrix :4x4 DH matrix
*/
    Matrix A;
    vector<double> row;

    for(int i=0; i<4; i++)
        row.push_back(0);

    for(int i=0; i<4; i++)
        A.push_back(row);

    A[0][0] = cos(theta);
    A[0][1] = -sin(theta) * cos(alpha);
    A[0][2] = sin(theta) * sin(alpha);
    A[0][3] = a * cos(theta);
    A[1][0] = sin(theta);
    A[1][1] = cos(theta) * cos(alpha);
    A[1][2] = -cos(theta) * sin(alpha);
    A[1][3] = a * sin(theta);
    A[2][1] = sin(alpha);
    A[2][2] = cos(alpha);
    A[2][3] = d;
    A[3][3] = 1;

    return A;
}

Matrix Jacobian(vector<double> l, vector<double> thetas){
/*
*   Inputus:
*       - vector<double> l: lenghts of 3 DoF RRR robot
*       - vector<double> thetas: angles of 3 DoF RRR robot
*
*   Output:
*       - Matrix : 3x3 Jacobian matrix
*/
    Matrix J;
    vector<double> row;

    for(int i=0; i<3; i++)
        row.push_back(0);

    for(int i=0; i<3; i++)
        J.push_back(row);

    J[0][0] = -l[2]*sin(thetas[0])*cos(thetas[1])*cos(thetas[2]) + l[2]*sin(thetas[0])*sin(thetas[1])*sin(thetas[2]) - l[1]*sin(thetas[0])*cos(thetas[1]);
    J[0][1] = -l[2]*cos(thetas[0])*sin(thetas[1])*cos(thetas[2]) - l[2]*cos(thetas[0])*cos(thetas[1])*sin(thetas[2]) - l[1]*cos(thetas[0])*sin(thetas[1]);
    J[0][2] = -l[2]*cos(thetas[0])*cos(thetas[1])*sin(thetas[2]) - l[2]*cos(thetas[0])*sin(thetas[1])*cos(thetas[2]);
    J[1][0] = l[2]*cos(thetas[0])*cos(thetas[1])*cos(thetas[2]) - l[2]*cos(thetas[0])*sin(thetas[1])*sin(thetas[2]) + l[1]*cos(thetas[0])*cos(thetas[1]);
    J[1][1] = -l[2]*sin(thetas[0])*sin(thetas[1])*cos(thetas[2]) - l[2]*sin(thetas[0])*cos(thetas[1])*sin(thetas[2]) - l[1]*sin(thetas[0])*sin(thetas[1]);
    J[1][2] = -l[2]*sin(thetas[0])*cos(thetas[1])*sin(thetas[2]) - l[2]*sin(thetas[0])*sin(thetas[1])*cos(thetas[2]);
    J[2][0] = 0;
    J[2][1] = l[1]*cos(thetas[1]) - l[2]*sin(thetas[1])*sin(thetas[2]) + l[2]*cos(thetas[1])*cos(thetas[2]);
    J[2][2] = l[2]*cos(thetas[1])*cos(thetas[2]) - l[2]*sin(thetas[1])*sin(thetas[2]);

    return J;
}

Matrix IK(Matrix &Path, vector<double> &lenghts){
/*
*   Inputs:
*       - Matrix Path: each row represents one point to be visited
*       - vector<double> lenghts: lenght of each link in the RRR robot
*
*   Outputs:
*       - Matrix : set of angles to perform travel
*/
    double epsilon = 1e-12;          // Set an error term
    int maxIter = 1000;             // Max of iterations
    vector<double> current_angle;   // vector with current angles (initial guess)
    vector<double> new_angle;      // vector with past calculated angles
    vector<double> current_point;   // vector with current point
    Matrix angles;                  // Matrix with calculated angles
    Matrix J, J_inv;
    vector<double> xd;

    cout << "Visited points:" << endl;
    for(int p=0; p<Path.size(); p++){
        // Set the initial angles (initial guess)
        if(p==0){
            for(int q=0; q<3; q++)
                current_angle.push_back(PI/2);
                angles.push_back(current_angle);
            }
        else
            current_angle = angles.back();

        xd = Path[p];               // Set the disired point
        for(int i=0; i<maxIter; i++){
            Matrix A1, A2, A3, T0_2, T0_3;

            A1 = DH(0,PI/2.0,lenghts[0],current_angle[0]);
            A2 = DH(lenghts[1],0.0,0.0,current_angle[1]);
            A3 = DH(lenghts[2],0.0,0.0,current_angle[2]);

            T0_2 = AXB(A1,A2);
            T0_3 = AXB(T0_2,A3);

            current_point.clear();
            for(int q=0; q<3; q++)
                current_point.push_back(T0_3[q][3]);

            vector<double> e = vectSubs(xd,current_point);

            double norme = normV(e);
            if(norme < epsilon){
                break;
            }
            else{
                J = Jacobian(lenghts,current_angle);
                J_inv = Inverse(J,3,3);

                vector<double> Je = matXvect(J_inv,e);
                new_angle = vectAdd(current_angle,Je);
                for(int q=0; q<3; q++)
                    new_angle[q] = fmod(new_angle[q],2*PI);
                angles.push_back(new_angle);
                current_angle = new_angle;
            }
        }

        if(p%10 == 0){
            for(int q=0; q<3; q++)
                cout << current_point[q] << "  ";
            cout << endl;
        }
    }
    cout << endl;
    return angles;
}

vector<double> matXvect(Matrix A,vector<double> V0){
/*
*   Inputs:
*       - Matrix A: Matrix A
*       - vector<double> V0: vector V0
*
*    Outputs:
*       - vector<double> : multiplication A*V0
*/

    vector<double> V1;

    for(int i=0; i<A.size(); i++){
        double sum = 0;
        for(int j=0; j<V0.size(); j++)
            sum += A[i][j]*V0[j];
        V1.push_back(sum);
    }

    return V1;
}

Matrix AXB(Matrix A,Matrix B){
/*
*   Inputs:
*       - Matrix A: Matrix A
*       - Matrix B: Matrix B
*
*    Outputs:
*       - Matrix : multiplication A*B
*/
    Matrix C;
    if(A[0].size() == B.size()){
        vector<double> row;

        for(int i=0; i<A.size(); i++){
            for(int j=0; j<B[0].size(); j++){
                double val = 0;
                for(int k=0; k<A[0].size(); k++)
                    val += A[i][k] * B[k][j];
                row.push_back(val);
            }
            C.push_back(row);
            row.clear();
        }
    }
    else{
        cout << "Error: A and B are NOT compatible\nMultiplication can not be done\n";
    }

    return C;
}


vector<double> vectSubs(vector<double> a, vector<double> b){
/*
*   Inputs:
*       - vector<double> a: vector a
*       - vector<double> a: vector b
*
*    Outputs:
*       - vector<double> : sum a - b
*/
    vector<double> c;
    double val;
    for(int i=0; i<a.size(); i++){
        val = a[i] - b[i];
        c.push_back(val);
    }

    return c;
}

vector<double> vectAdd(vector<double> a, vector<double> b){
/*
*   Inputs:
*       - vector<double> a: vector a
*       - vector<double> a: vector b
*
*    Outputs:
*       - vector<double> : sum a + b
*/
    vector<double> c;
    double val;
    for(int i=0; i<a.size(); i++){
        val = a[i] + b[i];
        c.push_back(val);
    }

    return c;
}

double normV(vector<double> V){
/*
*   Intputs:
*       - vector<double> V: vector to calculate its norm
*       - int n: number of entries in vector V
*
*   Outputs:
*       - double : euclidean norm of vector V
*/
    double norm = 0;
    for(int i=0; i<V.size(); i++) norm += V[i]*V[i];

    norm = sqrt(norm);

    return norm;
}

Matrix Inverse(Matrix A, int n, int m){
/*
*   Inputs:
*       - Matrix A: Pointer to the matrix A
*       - int n: number of rows in A
*       - int m: number of columns in A
*
*   Outputs:
*       - Matrix A_inv: Inverse matrix of A
*/
    Matrix A_inv;
    vector<double> row;
    for(int i=0; i<m; i++)
        row.push_back(0);

    for(int i=0; i<n; i++)
        A_inv.push_back(row);

    Matrix L;
    for(int i=0; i<=m; i++)     // Matrrix L has an extra column so it could be transformed to an augmented matrix
        row.push_back(0);

    for(int i=0; i<n; i++)
        L.push_back(row);

    Matrix U;
    for(int i=0; i<=m; i++)     // Matrrix L has an extra column so it could be transformed to an augmented matrix
        row.push_back(0);

    for(int i=0; i<n; i++)
        U.push_back(row);

    // A_inv is first set to the identity matrix
    for(int i=0; i<n;i++){
        for(int j=0; j<n; j++){
            if(i == j) A_inv[i][j] = 1.0;
            else A_inv[i][j] = 0.0;
        }
    }

    // First evaluates if this process can be done, A has to be a squared matrix
    if(n == m){
        // Get LU factorization from A
        // Initialize 1's and 0's in L and U
        for(int i=0;i<n;i++){
            for(int j=0;j<n;j++){
                if(i<j) L[i][j] = 0.0;
                else if(i>j) U[i][j] = 0.0;
                else L[i][j] = 1.0;
            }
        }

        int swap_occured = 0;   // If a swap between rows is nedded this counter tells how many have accurred

        for(int j=0;j<n;j += 1+swap_occured){
            double sum = 0;
            swap_occured = 0;

            // First row of U and first column of L are calculated
            if(j==0){
                for(int j=0;j<n;j++) U[0][j] = A[0][j];
                for(int i=1;i<n;i++) L[i][0] = A[i][0]/U[0][0];
            }

            // Use of Doolittle Algorithm
            for(int i=1; i<=j;i++){
                for(int k=0;k<=i-1;k++){
                    sum += L[i][k]*U[k][j];
                }
                U[i][j] = A[i][j]-sum;
                sum = 0;
            }

            for(int i=j+1;i<n;i++){
                for(int k=0;k<=j-1;k++){
                    sum += L[i][k]*U[k][j];
                }
                L[i][j] = (A[i][j]-sum)/U[j][j];
                sum = 0;
            }
            //////////////////////////////

            // Swap between rows if element on U[j][j]==0
            if(U[j][j] == 0){
                swap_occured--;

                // Swap of rows in matrix A
                double row_aux[n+1];
                for(int i=0;i<n;i++){
                    row_aux[i] = A[j+1][i];
                    A[j+1][i] = A[j][i];
                    A[j][i] = row_aux[i];
                }

                // Swap of rows in matrix U
                for(int i=0;i<n;i++){
                    row_aux[i] = U[j+1][i];
                    U[j+1][i] = U[j][i];
                    U[j][i] = row_aux[i];
                }

                // Swap of rows in matrix L
                for(int i=0;i<n;i++){
                    row_aux[i] = L[j+1][i];
                    L[j+1][i] = L[j][i];
                    L[j][i] = row_aux[i];
                }

                // Swap of rows in A_inv
                for(int i=0;i<n;i++){
                    row_aux[i] = A_inv[j+1][i];
                    A_inv[j+1][i] = A_inv[j][i];
                    A_inv[j][i] = row_aux[i];
                }

                // 0's and 1's are restored
                for(int i=0;i<n;i++){
                    for(int j=0;j<n;j++){
                        if(i<j) L[i][j] = 0.0;
                        else if(i>j) U[i][j] = 0.0;
                        else L[i][j] = 1.0;
                    }
                }
            }
        }

        //For each column in A_inv = I
        vector<double> x;       // Array where solution to x will be stored LUx = b
        vector<double> y;          // Array where solution to y will be stored Ly = b
        for(int i=0; i<n; i++){
            x.push_back(0);
            y.push_back(0);
        }

        for(int j=0;j<n;j++){

            for(int i = 0; i<n; i++) L[i][n] = A_inv[i][j];

            solver_TInfAug(L,n,y);

            // Take y as solution to augmented matrix Ux=y and solve for x
            for(int i = 0; i<n; i++) U[i][n] = y[i];

            solver_TSupAug(U,n,x);

            // Replace x in column j of A_inverse
            for(int i =0;i<n;i++) A_inv[i][j] = x[i];
        }
    }

   return A_inv;
}

void solver_TSupAug  (Matrix &A, int n,vector<double> &x){
/*
*   Inputs:
*       - Matrix A: pointer to augmented matrix
*       - int n: mumber of rows of matrix A
*
*   Output:
*       - vector<double> x: x solutions to augmented lower matrix A
*/
   double sum = 0;

    x[n-1] = A[n-1][n]/A[n-1][n-1];

    for(int i=n-2;i>=0;i--){
            sum = 0;
            for(int j=n-1;j>i;j--) sum += A[i][j]*x[j];
            x[i] = (A[i][n] - sum)/A[i][i];
    }
}

void solver_TInfAug(Matrix &A, int n,vector<double> &x){
/*
*   Inputs:
*       - Matrix A: pointer to augmented matrix
*       - int n: mumber of rows of matrix A
*
*   Output:
*       - vector<double> x: x solutions to augmented lower matrix A
*/
    double sum = 0;

    x[0] = A[0][n]/A[0][0];

    for(int i=1;i<n;i++){
            sum = 0;
            for(int j=0;j<i;j++) sum += A[i][j]*x[j];
            x[i] = (A[i][n] - sum)/A[i][i];
    }
}
