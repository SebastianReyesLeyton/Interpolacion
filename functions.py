"""
Name File: Functions.py 
Version: 1.0
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
    matmul     # Realize the matricial multiplication
)
import matplotlib.pyplot as plt

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

    return array(ans)

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
    Input:
    Description:
    Output:
    """
    ans, i, j, flag = None, col, col, M[col][col]
    n, m = len(M), len(M[0])

    while ( not flag and i < n ):
        j = col
        while ( not flag and j < m ):
            flag = M[i][j]
            if not flag: j += 1

        if not flag: i += 1
    ans = (i, j)
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
    #print(M[col][col])
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

def upperTriangular(A, b):
    """
    Input: The matrix A and the vector result b.
    Description: This function implement the Successive substitution backwards. Which is used 
                 when the matrix has The upper Triangular form.
    Output: The x values that resolve the system.
    """
    print(b)
    print(A)
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


def mult(p1, p2):
    ans = []
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
    for i in range(n):
        ans += p[i][0]*(t**i)
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

def graphics(polynomial, data, name):

    ans = proof(data[0], polynomial, len(polynomial) - 1)    
    
    plt.figure(name)
    plt.title(name)
    plt.plot(data[0], ans)
    plt.plot(data[0], data[1])
    plt.show()