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
    lowerTriangular,
    eval_polynomial,
    mult,
    sumPolynomials,
    solve,
    graphics
)

def polynomial(data):
    """
    Input: The matrix data of 2 x n data, which contain the t values on the first row and y values 
           in the last one.
    Description: Implementation of polynomial interpolation methods.
    Output: The solve to the dataset data.
    """
    ans, n = [], len(data[0])

    A = adjustment_matrix(n-1, n, data[0])
    b = array([data[1]])
    b = b.T

    print(A)
    print(b)
    #ans = gaussianElimination(A, b)
    ans = solve(A, b)

    return ans

def lagrange(data):
    """
    Input: The matrix data of 2 x n data, which contain the t values on the first row and y values 
           in the last one.
    Description: This function calls lagrange function, which calulate l_j(t) polynomials, and after multiply that polynomials
                 by its own y value, in others words, this function create the polynomial of the form:

                        p_(n-1)(t) = y_1 * l_1(t) + y_2 * l_2(t) + ... + y_(n-1) * l_(n-1)(t)

                Then of calculate all values, this function combine all polynomials and return the p_(n-1) polynomial,
                where n is the size of dataset.
    Output: The polynomial that is obtained by use lagrange method.
    """
    ans, n = [], len(data[0])
    ans = [ [ u*data[1][i] for u in lagrangeFunction(i, data[0])] for i in range(n) ]

    pt = [ 0 for _ in range(n) ]
    for i in range(n):
        for j in range(n):
            pt[i] += ans[j][i]

    ans = pt
    
    return array([ans]).T

def newton(data):
    """
    Input: The matrix data of 2 x n data, which contain the t values on the first row and y values 
           in the last one.
    Description:
    Output: The polynomial that is obtained by use newton method.
    """

    ans, n = array([ [ 0 for _ in range(len(data[0])) ] for _ in range(len(data[0])) ]), len(data[0])
    p = []

    for i in range(n): ans[i][0] = 1

    for i in range(n):

        p = newtonFunction(i, data[0])
        for j in range(len(p)):
            tmp = eval_polynomial(p[j], data[0][i], len(p[j])-1)
            ans[i][j+1] = tmp              
    
    ans = lowerTriangular(array(ans), array([data[1]]).transpose() )

    tmp = [ans[0][0]]
    for i in range(len(p)):
        p[i] = mult(p[i], [ans[i+1][0]])
        tmp = sumPolynomials(tmp, p[i])

    ans = tmp

    return array([ans]).T

def piecewiseLinear(data):
    """
    Input: The matrix data of 2 x n data, which contain the t values on the first row and y values 
           in the last one.
    Description:
    Output: The matrix of parameters for dataset.
    """

    ans = []

    for i in range(len(data[0])-1):

        p1 = [data[0][i], data[0][i+1]]
        A = adjustment_matrix(1, 2, p1)
        b = array([[data[1][i]],[data[1][i+1]]])
         
        ans.append(list(map(float,solve(A, b))))

    return array(ans)