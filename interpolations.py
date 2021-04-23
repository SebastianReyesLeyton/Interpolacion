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




