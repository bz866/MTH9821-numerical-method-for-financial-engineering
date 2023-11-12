// This is a program written as part of Baruch MFE course MTH 9821.
// 
// Author: Kevin Griffin
// 
// Revision: A
// Release Date: 11/10/2023
// Comments: Initial release.
// 
// Program Objective: NLA homework #9.
// 


#include "Decomposition.hpp"
#include "jacobi_solver.hpp"
#include <iostream>
#include <Dense>

void test_decomposition()
{
    // This function will test forward_subst, backward_subst, lu_row_pivoting, lu_no_pivoting, and cholesky. All test matrices
    // are taken from Dan's handouts to verify the solutions.
    //

    // Define test matrices for LU and Cholesky decomposition.
    mat A(4, 4);
    mat A2(5, 5);
    mat A3(4, 4);
    
    // Define a known solution x.
    vec x_true(4);
    
    // Perform LU decomposition without pivoting.
    A << 2, -1, 3, 0,
        -4, 5, -7, -2,
        -2, 10, -4, -7,
        4, -14, 8, 10;
    
    mat L, U;
    std::tie(L, U) = lu_no_pivoting(A);
    std::cout << "L: " << std::endl << L << std::endl;
    std::cout << "U: " << std::endl << U << std::endl;

    std::cout << std::endl;

    // Perform LU decomposition with row pivoting.
    A2 << 1, 2, -7, -1.5, 2,
        4, 4, 0, -6, -2,
        -2, -1, 2, 6, 2.5,
        0, -2, 2, -4, 1,
        2, 1, 13, -5, 3.5;

    permutation P;
    mat L_piv, U_piv;
    std::tie(P, L_piv, U_piv) = lu_row_pivoting(A2);
    //std::cout << "P: " << std::endl << P << std::endl;
    std::cout << "L_piv: " << std::endl << L_piv << std::endl;
    std::cout << "U_piv: " << std::endl << U_piv << std::endl;

    std::cout << std::endl;

    // Perform Cholesky decomposition.
    A3 << 9, -3, 6, -3,
        -3, 5, -4, 7,
        6, -4, 21, 3,
        -3, 7, 3, 15;

    mat U_chol = cholesky(A3);
    std::cout << "U_chol: " << std::endl << U_chol << std::endl;

    std::cout << std::endl;

    // Define a lower triangular matrix L (3x3 for example).
    mat L_sub(3, 3);
    L_sub << 1, 0, 0,  // Lower triangular matrix.
          2, 1, 0,
          3, 4, 1;

    // Define an upper triangular matrix U (3x3 for example).
    mat U_sub(3, 3);
    U_sub << 1, 2, 3,  // Upper triangular matrix.
          0, 1, 4,
          0, 0, 1;

    // Define a right-hand side vector b.
    vec b(3);
    b << 5,  // Vector b
          6,
          7;

    // Perform forward substitution with L and b to solve for y.
    vec y = forward_subst(L_sub, b); // Expected: [5, -4, 8] per Octave.

    // Perform backward substitution with U and y to solve for x.
    vec x = backward_subst(U_sub, b); // Expected: [28, -22, 7] per Octave.

    // Verify the solutions.
    std::cout << "y (solution of Ly = b):" << std::endl << y << std::endl;
    std::cout << "x (solution of Ux = b):" << std::endl << x << std::endl;
}


void test_jacobi_solver()
{
    // Define a matrix A and vector b representing the linear system Ax = b
    mat A(3, 3);
    A << 4, 1, 1,
         1, 3, -1,
         1, -1, 3;
    vec b(3);
    b << 6, 8, 7;

    // Initial guess x_0 (can be zeros or any other guess)
    vec x_0 = vec::Zero(3);

    // Tolerance for stopping criterion
    double tolerance = 0.001;

    // Solve the system using the Jacobi method
    vec x;
    unsigned int iterations;
    std::tie(x, iterations) = jacobi(A, b, x_0, tolerance, consecutive);

    // Output the solution
    std::cout << "Jacobi with Consecutive Stopping Criteria" << std::endl;
    std::cout << "Solution vector x:\n" << x << std::endl;
    std::cout << "Iterations: " << iterations << std::endl;

    std::tie(x, iterations) = jacobi(A, b, x_0, tolerance, residual);
    std::cout << "Jacobi with Residual Stopping Criteria" << std::endl;
    std::cout << "Solution vector x:\n" << x << std::endl;
    std::cout << "Iterations: " << iterations << std::endl;
}

int main()
{
    test_decomposition();
    test_jacobi_solver();
    return 0;
}