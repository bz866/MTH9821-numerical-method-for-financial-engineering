// This is a program written as part of Baruch MFE course MTH 9821.
// 
// Author: Kevin Griffin
// 
// Revision: A
// Release Date: 11/10/2023
// Comments: Initial release.
// 
// Program Objective: Implements functions needed for performing LU and Cholesky decompositions using the Eigen library.
// 

#include "Decomposition.hpp"

// Substitution functions. //
vec forward_subst(const mat& L, const vec& b)
{
    // Forward substitution algorithm for use with LU solvers.
	// Inputs:
	//	* L: Non-singular, lower triangular matrix of size n.
    //  * b: Column vector of size n (RHS).
	// Outputs:
	//	* x: Solution to Lx = b.
	//
	// Reference: Stefanica pg. 40.
    
    long n = b.size();  // Size of RHS vector.
    vec x = vec(n);     // Initialize solution vector as zeros.

 
    x(0) = b(0) / L(0, 0);
    for (long j = 1; j < n; ++j)
    {
        double sum = 0.0;
        for (long k = 0; k < j; ++k)
        {
            sum += L(j, k) * x(k);
        }

        x(j) = (b(j) - sum) / L(j, j);
    }

    return x;
}

vec backward_subst(const mat& U, const vec& b)
{
    // Backward substitution algorithm for use with LU solvers.
	// Inputs:
	//	* U: Non-singular, upper triangular matrix of size n.
    //  * b: Column vector of size n (RHS).
	// Outputs:
	//	* x: Solution to Ux = b.
	//
	// Reference: Stefanica pg. 45.
    
    long n = b.size();  // Size of RHS vector.
    vec x = vec(n);     // Initialize solution vector as zeros.

    // Start with the last element and move upwards.
    x(n - 1) = b(n - 1) / U(n - 1, n - 1);
    
    for (long j = n - 2; j >= 0; --j)
    {
        double sum = 0;
        for (long k = j + 1; k < n; ++k)
        {
            sum += U(j, k) * x(k);
        }
        x(j) = (b(j) - sum) / U(j, j);
    }

    return x;
}

// LU Decomposition functions. //
std::tuple<mat, mat> lu_no_pivoting(mat A)
{
	// Computes the LU decomposition of a matrix without row-pivoting. Conditions for input matrix:
	//	* Non-singular.
	//	* All leading principle minors are non-zero.
	//	* Full rank.
	//	* Typical examples: SPD matrices, diagonally dominant matrices.
	// Inputs:
	//	* A: Matrix of size n x n.
	// Outputs:
	//	* L: Lower triangular matrix with 1's on the main diagonal.
	//	* U: Upper triangular matrix such that A = LU.
	//
	// Reference: Stefanica pg. 54.

	long n = A.rows();				// Size of square matrices.

    mat L = mat::Identity(n, n);	// Lower triangular matrix -- initialize with 1's on diagonal.
    mat U = mat::Zero(n, n);		// Upper triangular matrix.

	for (long i = 0; i < n - 1; ++i)
	{
        for (long k = i; k < n; ++k)
		{
            // Compute row i of U.
            U(i, k) = A(i, k);
            
			// Compute column i of L.
            if (U(i, i) != 0) // Check to prevent division by zero.
			{
                L(k, i) = A(k, i) / U(i, i);
            }
        }
        
        for (long j = i + 1; j < n; ++j)
		{
            for (long k = i + 1; k < n; ++k)
			{
                // Modify the future elements to be decomposed.
                A(j, k) -= L(j, i) * U(i, k);
            }
        }
    }

    // Set the bottom-right element of L to 1 and U to the corresponding element of A.
    L(n - 1, n - 1) = 1;
    U(n - 1, n - 1) = A(n - 1, n - 1);

    return std::make_tuple(L, U);
}

std::tuple<permutation, mat, mat> lu_row_pivoting(mat A)
{
	// Computes the LU decomposition of a matrix with row-pivoting. Conditions for input matrix:
	//	* Non-singular.
	//	* Square (n x n).
	// Inputs:
	//	* A: Matrix of size n x n.
	// Outputs:
	//	* P: Permutation matrix stored as a vector of diagonal entries.
	//	* L: Lower triangular matrix with 1's on the main diagonal.
	//	* U: Upper triangular matrix such that A = LU.
    //
	// Reference: Stefanica pg. 70.

	long n = A.rows();              // Size of square matrices.
    permutation P(n);               // Pemutation matrix -- initialize with zeros.
    mat L = mat::Identity(n, n);    // Lower triangular matrix -- initialize as identity.
    mat U = mat::Zero(n, n);		// Upper triangular matrix.

    // Initialize the permutation matrix P to the identity permutation.
    P.setIdentity();

    for (long i = 0; i < n; ++i)
    {
        // Finding the pivot.
        long imax = i;
        double maxA = -10e8;
        for (long k = i; k < n; ++k)
        {
            if (std::abs(A(k, i)) > maxA)
            {
                maxA = std::abs(A(k, i));
                imax = k;
            }
        }
        
        // Row interchange for the pivot.
        if (imax != i)
        {
            A.row(i).swap(A.row(imax));
            P.applyTranspositionOnTheRight(i, imax);
            if (i > 0)
            {
                L.block(i, 0, 1, i).swap(L.block(imax, 0, 1, i));
            }
        }

        for (long j = i; j < n; ++j)
        {
            L(j, i) = A(j, i) / A(i, i);
            U(i, j) = A(i, j);
        }            
        for (long j = i + 1; j < n; ++j)
        {
            for (long k = i + 1; k < n; ++k)
            {
                A(j, k) -= L(j, i) * U(i, k);
            }
        }
    }

    // Finalizing L and U for last row.
    L(n - 1, n - 1) = 1;                // Ensure the last diagonal element of L is 1.
    U(n - 1, n - 1) = A(n - 1, n - 1);  // The last element of U is the bottom-right element of A.

    return std::make_tuple(P, L, U);
}

// Cholesky Decomposition function. //
mat cholesky(mat A)
{
    // Computes the Cholesky decomposition of a matrix. Conditions for input matrix:
	//	* Non-singular.
	//	* SPD.
	// Inputs:
	//	* A: Matrix of size n x n.
	// Outputs:
	//	* U: Upper triangular matrix such that A = U^T*U.
    //
	// Reference: Stefanica pg. 170.

    long n = A.rows();          // Size of square matrices.
    mat U = mat::Zero(n, n);    // Initialize U as a zero matrix of size n x n.

    for (long i = 0; i < n; ++i)
    {
        // Compute the diagonal element.
        U(i, i) = std::sqrt(A(i, i));
        
        for (long k = i + 1; k < n; ++k)
        {
            U(i, k) = A(i, k) / U(i, i);  // Compute row i of U.
        }

        for (long j = i + 1; j < n; ++j)
        {
            for (long k = j; k < n; ++k)
            {
                // Update the remaining elements.
                A(j, k) -= U(i, j) * U(i, k);
            }
        }
    }

    // Compute the last diagonal element of U.
    U(n - 1, n - 1) = std::sqrt(A(n - 1, n - 1));

    return U;
}