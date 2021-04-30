"""
Name File: Functions.py 
Version: 2.5
Description: A file with external functions used on interpolations.py 
             classe

Participants:
- Juan Sebastián Reyes
- Alexander Castro
"""

from numpy import (
    transpose, # Calculate the transpose of a matrix
    array,     # Type of data to allow use numpy functions
    dot,       # Calculate the point product
    identity,  # Generate the identity matrix
    matmul,    # Realize the matricial multiplication
    linspace,  # Return a quantity num of points between start and stop
    mean,      # Calculate the mean
    std        # Calculate de standard desviation
    
)
from numpy.linalg import solve  # Resolve the Ax = b linear system
import matplotlib.pyplot as plt

# Improvements

def binary_exponentiation(n, k):
    """
    Input: The base n and the exponent k
    Description: This function implement and return a strategy to calculate to
                 a number n elevate to some exponent k. This tecnique use the 
                 concept of calculate the half of k and after obtain the result
                 depend if the k is even or odd is performed different calculations.
                 Formally is:
                        if k = 0, return 1
                        if k = 1, return n
                        else
                            if k mod 2 = 0, is calculated n^(k/2) and this result is multiply by itself.
                            else, is calculated n^(k/2) and this result is multiply by itself and n.
    Output: The value of n to k
    """
    ans = 1
    if (k == 1): ans = n
    elif (k > 1):
        ans = binary_exponentiation(n, k // 2)
        if ( k % 2 ):  ans = ans*ans*n
        else: ans = ans*ans
    return ans

# Creation of the system

def adjustment_matrix(n, m, t):
    """
    Input: The degree of the polynomial n, the amount of data m, and 
           the values t.
    Description: This function create and return a adjustment matrix of 
                 n degree, which is defined as:
                                    A of m x n 
                 where A_ij = t_i^j , i and j between 0 .. n-1.
    Output: The A matrix.
    """
    ans = [ [ 0 for _ in range(n+1)] for _ in range(m) ]

    for i in range(m):
        for j in range(n+1):
            ans[i][j] = t[i]**j

    return array(ans, dtype='float')

# Matrix operations

def permutations(M, i, j, row):
    """
    Input: The matrix M, columns/rows i and j, and a flag row that specify whether the permutation is in rows or columuns.
    Description:
    Output:
    """
    ans = identity(len(M))
    ans[i][j] = ans[i][i]; ans[i][i] = 0
    ans[j][i] = ans[j][j]; ans[j][j] = 0
    if row: M = matmul(ans, M)
    else: M = matmul(M, ans)
    ans = M
    return ans

def searchPivot(M, col):
    """
    Input: The matrix M and the index col
    Description: This functions seach a position in the matrix that its value is not 0,
                 then return that position (i, j).
                 Note: This function search between col <= i, j < n 
    Output: The position i and j, where i represent the row position and j column position.
    """
    ans, i, j, flag = None, col, col, M[col][col]
    n, m = len(M), len(M[0])

    while ( not flag and i < n ):
        j = col
        while ( not flag and j < m ):
            flag = M[j][i]
            if not flag: j += 1

        if not flag: i += 1
    ans = (j, i)
    return ans

def eliminationMatrix(M, col):
    """
    Input: The matrix M and the column col
    Description: This functions implement calculate the elimination matrix by M in
                 the column col. In this process can occur to when we take a pivot this 
                 has value 0, so in this case, the permutation functions is used to calculate
                 a permutation of the matrix M, which the position M[col][col] is different to 0.
    Output: The elimination matrix to M and M with its changes.
    """
    
    ans = identity(len(M))
    if ( M[col][col] != 0 ):
        for i in range(col+1, len(M)):
            ans[i][col] = (-1 * M[i][col])/M[col][col]
    else:
        i, j = searchPivot(M, col)
        if (i > len(M) or j > len(M)):
            raise TypeError("The pivot can not be 0")
        else:
            if ( i == col and j != col ):
                M = permutations(M, col, j, 0)
            elif ( j == col and i != col ):
                M = array(permutations(M, col, i, 1))
            else:
                M = array(permutations(M, col, j, 0))
                M = array(permutations(M, col, i, 1))
            for i in range(col+1, len(M)):
                ans[i][col] = (-1 * M[i][col])/M[col][col]

    return array(ans), M

# Methods to resolve matrix systems

def upperTriangular(A, b):
    """
    Input: The matrix A and the vector result b.
    Description: This function implement the Successive substitution backwards. Which is used 
                 when the matrix has The upper Triangular form.
    Output: The x values that resolve the system.
    """
    n = len(b)
    
    ans = array([ [None] for _ in range(n) ])
    ans[n-1][0] = b[n-1][0] / A[n-1][n-1]

    for i in range(n-2, -1, -1):
        ans[i][0] = b[i][0]
        tmp = 0
        for j in range(i+1, n):
            tmp += A[i][j] * ans[j][0]

        ans[i][0] -= tmp
        ans[i][0] /= A[i][i]

    return ans

def lowerTriangular(A, b):
    """
    Input: The matrix A and the vector result b.
    Description: This function implement the Successive substitution forwards. Which is used when
                 the matrix has the Lower Triangular form.
    Output: The x values that resolve the system.
    """

    n = len(b)
    ans = [ [0] for _ in range(n) ]
    ans[0][0] = b[0][0] / A[0][0]

    for i in range(1, n):
        ans[i][0] = b[i][0]
        tmp = 0
        for j in range(i):
            tmp = tmp + A[i][j]*ans[j][0]

        ans[i][0] -= tmp
        ans[i][0] /= A[i][i]

    return ans

def gaussianElimination(M, b):
    """
    Input:
    Description:
    Output:
    """

    ans, Me, br = M, [], b
    
    for i in range(len(M) - 1):
        Me, ans = eliminationMatrix(ans, i)
        ans = matmul(Me, ans)
        br = matmul(Me, br)

    ans = upperTriangular(ans, br)
    return ans

# Polynomial operations

def mult(p1, p2):
    """
    Input: The polynomial p1 and p2
    Description: This function perform the multiplication between p1 and p2
    Output: The result of multiply p1 by p2
    """
    ans = [ 0 for _ in range(len(p1)-1 + len(p2)-1 + 1) ]

    for i in range(len(p1)):
        for j in range(len(p2)):
            ans[i+j] += p1[i]*p2[j]

    return ans

def sumPolynomials(p1, p2):
    """
    Input: The polynomials p1 and p2.
    Description: This function sum both polynomials.
    Output: The result of sum p1 and p2.
    """
    ans = []

    if (len(p1) > len(p2)):
        n, m = len(p1), len(p2)
        for i in range(n):
            if (i < m): ans.append(p1[i] + p2[i])
            else: ans.append(p1[i])
    else:
        n, m = len(p2), len(p1)
        for i in range(n):
            if (i < m): ans.append(p1[i] + p2[i])
            else: ans.append(p2[i])

    return ans

# Auxiliar function by lagrange and newton interpolation

def lagrangeFunction(j, t):
    """
    Input: The index j and the array with t values.
    Description: This function calculate the polynomial create by l_j(t) of lagrange function, where
                 this is define as:
                                     __
                            l_j(t) = || (t - t_k)
                                    _____________
                                     __
                                     || (t_j - t_k)
                 
                 where both products start with k=1, but k is always different of j, and finish when k > n
    Output: The polynomial that is produced by lagrange funtion.
    """
    ans, tmp,  n = [1], 1, len(t)

    for k in range(n):

        if ( k != j ):

            # Calculate the numerator
            p = [-t[k], 1]
            ans = mult(p, ans)

            # Calculare the denominator
            tmp *= (t[j] - t[k])

    for i in range(len(ans)): ans[i] = ans[i] / tmp

    return ans

def newtonFunction(j, t):
    """
    Input: The index j and the array with the t values
    Description: This function calculate the Newton phi function, which is define as:
                               __
                    phi_j(t) = || (t - t_k)

                 Where this product start with k = 1 and finish when k > j-1. 
    Output: The polynomial obtained by the product 
    """
    ans, tmp = [ [] for _ in range(j) ], [1]

    for k in range(j):
        p = [-t[k], 1]
        tmp = mult(p, tmp)
        ans[k] = tmp
        
    return ans

# Graphics

def eval_polynomial(p, t, n):
    """
    Input: The polynomial p, the value of t variable and the degree n of the polynomial p.
    Description: This functions assign t value to t variable in the polynomial p of degree n,
                 which is described as follows:
                            p(t) = x_1 + x_2*t + ... + x_n*t^{n-1}
                 This result obtain when is evaluate t in the polynomial p is returned.
    Output: Result to assess t in the linear polynomial p.
    """
    ans = 0
    for i in range(n+1):
        ans += p[i]*(t**i)
    return ans

def proof(test, p, n):
    """
    Input: The array test which have the t values, the polynomial p, and the degree of the polynomial n
    Description: This function realice all evaluation of t values in p polynomial through eval_polynomial 
                 function and return those answers.
    Output: The evaluation of all t values in polynomial p.
    """
    ans = []
    for t in test:
        ans.append(eval_polynomial(p, t, n))

    return ans

def graphics(polynomial, train, test, name):
    """
    Input: The polynomial, the train and test dataset and the name of the graphic
    Description: This function graphics train and test points, and the polynomial in the same
                 graphic, then show that graphic with the name pass by parameters.
    Output: None
    """ 
    print(len(polynomial))
    t = linspace(0, train[0][-1], 1000)
    ans = proof(t, polynomial, len(polynomial) - 1)    
    
    plt.figure(name)
    plt.title(name)
    plt.plot(t, ans)
    plt.plot(test[0], test[1], 'o')
    plt.plot(train[0], train[1], 'o')
    plt.legend(
        loc = 'upper right',
        labels=['Training', 'Polynomial','Test']
    )
    plt.xlabel(f'Quantity: {len(train[0])}')
    plt.grid()
    plt.show()

# Error function

def error(real, estimate):
    """
    Input: The list of real values and the list of estimate values.
    Description: This function calculate the average mean and the standard desviation of
                 absolute error between real and estimate values.
    Output: The mean and standard desviation of the error.
    """
    ans = [ r - e for r, e in zip(real, estimate)]
    return mean(ans), std(ans)