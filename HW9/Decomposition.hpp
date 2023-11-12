// This is a program written as part of Baruch MFE course MTH 9821.
// 
// Author: Kevin Griffin
// 
// Revision: A
// Release Date: 11/10/2023
// Comments: Initial release.
// 
// Program Objective: Defines functions needed for performing LU and Cholesky decompositions using the Eigen library.
// 

#ifndef Decomposition_HPP  // Conditional compilation
#define Decomposition_HPP

#include <iostream>
#include <tuple>
#include <Dense>     // Eigen package -- various functions for working with matrices. Source: https://eigen.tuxfamily.org/index.php?title=Main_Page.

// Eigen typedefs. //
typedef Eigen::VectorXd vec;
typedef Eigen::MatrixXd mat;
typedef Eigen::PermutationMatrix<-1 , -1 , unsigned long> permutation;

// Substitution functions. //
vec forward_subst(const mat& L, const vec& b);
vec backward_subst(const mat& U, const vec& b);

// LU Decomposition functions. //
std::tuple<mat, mat> lu_no_pivoting(mat A);
std::tuple<permutation, mat , mat> lu_row_pivoting(mat A);

// Cholesky Decomposition function. //
mat cholesky(mat A);

#endif