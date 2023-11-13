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
#include "gauss_siedel_solver.hpp"
#include "sor_solver.hpp"
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


void test_jacobi()
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
    std::cout << "\n\n" << std::endl;

    std::tie(x, iterations) = jacobi(A, b, x_0, tolerance, residual);
    std::cout << "Jacobi with Residual Stopping Criteria" << std::endl;
    std::cout << "Solution vector x:\n" << x << std::endl;
    std::cout << "Iterations: " << iterations << std::endl;
    std::cout << "\n\n" << std::endl;

}

void test_jacobi_banded()
{
    // Define a tridiagonal matrix A (bandwidth = 1)
    int n = 5; // Size of the matrix
    mat A = mat::Zero(n, n);
    A.diagonal() << 4, 3, 5, 3, 4; // Main diagonal
    A.diagonal(-1) << -1, -1, -1, -1; // Sub-diagonal
    A.diagonal(1) << -1, -1, -1, -1; // Super-diagonal

    // Define vector b
    vec b(n);
    b << 3, 3, 4, 3, 3;

    // Tolerance for convergence
    double tolerance = 0.001;
    int bandwidth = 1; // Tridiagonal matrix has bandwidth of 1

    vec x;
    unsigned int iterations;
    std::tie(x, iterations) = jacobiBanded(A, b, bandwidth, tolerance, consecutive);

    // Output the solution
    std::cout << "Jacobi Banded with Consecutive Stopping Criteria" << std::endl;
    std::cout << "Solution vector x:\n" << x << std::endl;
    std::cout << "Iterations: " << iterations << std::endl;
    std::cout << "\n\n" << std::endl;

    std::tie(x, iterations) = jacobiBanded(A, b, bandwidth, tolerance, residual);
    std::cout << "Jacobi Banded with Residual Stopping Criteria" << std::endl;
    std::cout << "Solution vector x:\n" << x << std::endl;
    std::cout << "Iterations: " << iterations << std::endl;
    std::cout << "\n\n" << std::endl;
}

void test_gauss_siedel()
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
    std::tie(x, iterations) = gauss_seidel(A, b, x_0, tolerance, consecutive);

    // Output the solution
    std::cout << "Gauss-Seidel with Consecutive Stopping Criteria" << std::endl;
    std::cout << "Solution vector x:\n" << x << std::endl;
    std::cout << "Iterations: " << iterations << std::endl;
    std::cout << "\n\n" << std::endl;

    std::tie(x, iterations) = gauss_seidel(A, b, x_0, tolerance, residual);
    std::cout << "Gauss-Seidel with Residual Stopping Criteria" << std::endl;
    std::cout << "Solution vector x:\n" << x << std::endl;
    std::cout << "Iterations: " << iterations << std::endl;
    std::cout << "\n\n" << std::endl;
}

void test_gauss_siedel_banded()
{
    // Define a tridiagonal matrix A (bandwidth = 1)
    int n = 5; // Size of the matrix
    mat A = mat::Zero(n, n);
    A.diagonal() << 4, 3, 5, 3, 4; // Main diagonal
    A.diagonal(-1) << -1, -1, -1, -1; // Sub-diagonal
    A.diagonal(1) << -1, -1, -1, -1; // Super-diagonal

    // Define vector b
    vec b(n);
    b << 3, 3, 4, 3, 3;

    // Tolerance for convergence
    double tolerance = 0.001;
    int bandwidth = 1; // Tridiagonal matrix has bandwidth of 1

    vec x;
    unsigned int iterations;
    std::tie(x, iterations) = gauss_seidel_banded(A, b, bandwidth, tolerance, consecutive);

    // Output the solution
    std::cout << "Gauss-Seidel Banded with Consecutive Stopping Criteria" << std::endl;
    std::cout << "Solution vector x:\n" << x << std::endl;
    std::cout << "Iterations: " << iterations << std::endl;
    std::cout << "\n\n" << std::endl;

    std::tie(x, iterations) = gauss_seidel_banded(A, b, bandwidth, tolerance, residual);
    std::cout << "Gauss-Seidel Banded with Residual Stopping Criteria" << std::endl;
    std::cout << "Solution vector x:\n" << x << std::endl;
    std::cout << "Iterations: " << iterations << std::endl;
    std::cout << "\n\n" << std::endl;
}

void test_sor()
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

    // Relaxation factor (omega), often set between 1 and 2
    double omega = 1.5;

    vec x;
    unsigned int iterations;
    std::tie(x, iterations) = sor(A, b, omega, tolerance, consecutive);

    // Output the solution
    std::cout << "SOR with Consecutive Stopping Criteria" << std::endl;
    std::cout << "Solution vector x:\n" << x << std::endl;
    std::cout << "Iterations: " << iterations << std::endl;
    std::cout << "\n\n" << std::endl;

    std::tie(x, iterations) = sor(A, b, omega, tolerance, residual);
    std::cout << "SOR with Residual Stopping Criteria" << std::endl;
    std::cout << "Solution vector x:\n" << x << std::endl;
    std::cout << "Iterations: " << iterations << std::endl;
    std::cout << "\n\n" << std::endl;
}

void test_sor_banded()
{
    // Define a tridiagonal matrix A (bandwidth = 1)
    int n = 5; // Size of the matrix
    mat A = mat::Zero(n, n);
    A.diagonal() << 4, 3, 5, 3, 4; // Main diagonal
    A.diagonal(-1) << -1, -1, -1, -1; // Sub-diagonal
    A.diagonal(1) << -1, -1, -1, -1; // Super-diagonal

    // Define vector b
    vec b(n);
    b << 3, 3, 4, 3, 3;

    // Tolerance for convergence
    double tolerance = 0.001;
    int bandwidth = 1; // Tridiagonal matrix has bandwidth of 1

    // Relaxation factor (omega)
    double omega = 1.5;

    vec x;
    unsigned int iterations;
    std::tie(x, iterations) = sor_banded(A, b, bandwidth, omega, tolerance, consecutive);

    // Output the solution
    std::cout << "SOR Banded with Consecutive Stopping Criteria" << std::endl;
    std::cout << "Solution vector x:\n" << x << std::endl;
    std::cout << "Iterations: " << iterations << std::endl;
    std::cout << "\n\n" << std::endl;

    std::tie(x, iterations) = sor_banded(A, b, bandwidth, omega, tolerance, residual);
    std::cout << "SOR Banded with Residual Stopping Criteria" << std::endl;
    std::cout << "Solution vector x:\n" << x << std::endl;
    std::cout << "Iterations: " << iterations << std::endl;
    std::cout << "\n\n" << std::endl;
}


void test_solve()
{
    std::cout << "Q6 Solve test" << std::endl;
    int n = 14; // Size of the matrix
    mat A = mat::Zero(n, n);
    A.diagonal() << 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2; // Main diagonal
    A.diagonal(-1) << -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1; // Sub-diagonal
    A.diagonal(1) << -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1; // Super-diagonal

    // Define vector b
    vec b(n);
    b << 0, 1, 4, 9, 16, 25, 36, 49, 64, 81, 100, 121, 144, 169;
    
    // Tolerance for convergence
    double tolerance = 0.000001;
    int bandwidth = 1; // Tridiagonal matrix has bandwidth of 1

    // Initial guess x_0
    vec x_0(n);
    x_0 << 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1;

    vec x;
    unsigned int iterations;

    // Q6.1 Jacobi iteration
    std::cout << "Q6.1 Solve the system using the Jacobi method" << std::endl;
    // Output the solution
    std::tie(x, iterations) = jacobiBanded(A, b, bandwidth, tolerance, residual);
    std::cout << "Jacobi Banded with Residual Stopping Criteria" << std::endl;
    std::cout << "Solution vector x:\n" << x << std::endl;
    std::cout << "Iterations: " << iterations << std::endl;
    std::cout << "\n\n" << std::endl;

    std::tie(x, iterations) = jacobiBanded(A, b, bandwidth, tolerance, consecutive);
    std::cout << "Jacobi Banded with Consecutive Stopping Criteria" << std::endl;
    std::cout << "Solution vector x:\n" << x << std::endl;
    std::cout << "Iterations: " << iterations << std::endl;
    std::cout << "\n\n" << std::endl;

    // Q6.2 GS iteration
    std::cout << "Q6.2 Solve the system using the GS method" << std::endl;
    // Output the solution
    std::tie(x, iterations) = gauss_seidel_banded(A, b, bandwidth, tolerance, residual);
    std::cout << "Gauss-Seidel Banded with Residual Stopping Criteria" << std::endl;
    std::cout << "Solution vector x:\n" << x << std::endl;
    std::cout << "Iterations: " << iterations << std::endl;
    std::cout << "\n\n" << std::endl;

    std::tie(x, iterations) = gauss_seidel_banded(A, b, bandwidth, tolerance, consecutive);
    std::cout << "Gauss-Seidel Banded with Consecutive Stopping Criteria" << std::endl;
    std::cout << "Solution vector x:\n" << x << std::endl;
    std::cout << "Iterations: " << iterations << std::endl;
    std::cout << "\n\n" << std::endl;

    // Q6.3 SOR iteration
    std::cout << "Q6.3 Solve the system using the SOR method" << std::endl;
    // Relaxation factor (omega)
    double omega = 1.15;
    // Output the solution
    std::tie(x, iterations) = sor_banded(A, b, bandwidth, omega, tolerance, residual);
    std::cout << "SOR Banded with Residual Stopping Criteria and omega = 1.15" << std::endl;
    std::cout << "Solution vector x:\n" << x << std::endl;
    std::cout << "Iterations: " << iterations << std::endl;
    std::cout << "\n\n" << std::endl;

    std::tie(x, iterations) = sor_banded(A, b, bandwidth, omega, tolerance, consecutive);
    std::cout << "SOR Banded with Consecutive Stopping Criteria and omega = 1.15" << std::endl;
    std::cout << "Solution vector x:\n" << x << std::endl;
    std::cout << "Iterations: " << iterations << std::endl;
    std::cout << "\n\n" << std::endl;

    // Q6.4 Change of omega on SOR
    std::cout << "Q6.4 The SOR method results for different omega's" << std::endl;
    for (float i = 1.02; i <= 1.98; i = i + 0.02)
    {
        omega = i;
        std::tie(x, iterations) = sor_banded(A, b, bandwidth, omega, tolerance, residual);
        std::cout << "omega = " << i << ", number to iterations = " << iterations << std::endl;
    }

}


int main()
{
    // Q1 Q2 Q3
    test_decomposition();
    // Q4
    test_jacobi();
    test_gauss_siedel();
    test_sor();
    // Q4 banded matrices
    test_jacobi_banded();
    test_gauss_siedel_banded();
    test_sor_banded();

    // Q6 
    test_solve();
    // 1. The residual-based and consecutive approximation stopping conditions result in identical solution and quite similar
    // convergence speed.
    // 2. Compare the convergence speed of the Jacobi method, the GS method and the SOR method:
    // SOR with a good choice of omega > GS > Jacobi > SOR with a bad choice of omega
    // 3. For the SOR method, when omega is near to 2, the convergence speed is rather low compared to the omega close to 1
    // 4. For the SOR method, when omega is near to 1, the solution and the convergence speed is similar to GS method (by definition).
    
    return 0;
}