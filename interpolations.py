"""
Name File: Interpolation.py
Version: 1.0
Description: Implementation of interpolation methods
             Polynomial, Lagrage, Newton and piecewise
Participants:
- Alexander Castro
- Juan Sebatián Reyes

"""

from functions import (
    adjustment_matrix,
    array,
    gaussianElimination,
    lagrangeFunction,
    newtonFunction,
    lowerTriangular
)

def polynomial(data):
    """
    Input: Matrix 2*n with data set
    Description: Implementation of polynomial interpolation methods
    Output: 
    """
    ans, n = [], len(data[0])

    A = adjustment_matrix(n-1, n, data[0])
    print(A)
    b = array([data[1]])
    b = b.T
    print(b)

    ans = gaussianElimination(A, b)
    return ans

def lagrange(data):
    """
    Input:
    Description:
    Output:
    """
    ans, n = [], len(data[0])
    ans = [ [ u*data[1][i] for u in lagrangeFunction(i, data[0])] for i in range(n) ]

    pt = [ 0 for _ in range(n-1) ]
    for i in range(n):
        for j in range(n):
            pt[i] += ans[j][i]

    ans = pt
    
    return ans

def newton(data):
    """
    Input:
    Description:
    Output:
    """

    ans, n = [ [ 0 for _ in range(len(data[0])) ] for _ in range(len(data[0])) ], len(data[0])

    for i in range(n):

        p = newtonFunction(i, data[0])
        for j in range(len(p)):
            ans[i][j] = p[j]

    ans = lowerTriangular(array(ans), array([data[1]]).transpose() )
