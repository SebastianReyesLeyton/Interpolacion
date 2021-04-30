"""
Name: main DataSet1.py
Version: 2.2
Participants:
    - Juan Sebastián Reyes
    - Alexander Castro
"""
from numpy import mean, std, float64
from pandas import read_csv, DataFrame
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from interpolations import (
    polynomial,
    lagrange,
    newton,
    piecewiseLinear
)
from functions import (
    graphics,
    error,
    proof
)
from time import time

def main():

    # Upload data of .csv file
    data = read_csv("./DataSet1.csv", header=None) # The time is expressed on minutes.

    # Obtain the array of time and its linked values
    T = data.drop(1, axis=1)
    y = data.drop(0, axis=1)

    # Obtain the first value of time
    d = T[0][0]

    # Transform the time array
    for i in range(T.shape[0]): T[0][i] -= d

    # The best jump 20
    jump = 2

    # Create the subset of train and test, and their corresponds outputs
    T_train, T_test, y_train, y_test = [T[0][0]] , [], [y[1][0]], []

    for i in range(1, T.shape[0]-1):

        if (i % (jump)):
            T_test.append(T[0][i])
            y_test.append(y[1][i])
        else:
            T_train.append(T[0][i])
            y_train.append(y[1][i])

    T_train.append(T[0][T.shape[0]-1])
    y_train.append(y[1][y.shape[0]-1])

    # Create the train matrix
    train = [list(map(int, T_train)), list(map(int, y_train))]
    
    # Create the test matrix
    test = [T_test, y_test]

    Data = [[-1, -0.5, 0.0, 0.5, 1], [1, 0.5, 0.0, 0.5, 2.0]]

    # Obtain the polynomial through normalize equations
    start = time()
    polynom = piecewiseLinear(train)
    stop = time()
    
    # Show the polynomial
    print(polynom)

    # Show the runtime of the interpolation method
    print(f"Time: {stop - start}")

    # If the piecewise method?
    if (len(polynom[0]) == 1):
        polynom = list(polynom.T[0])
        p = proof(test[0] , polynom, len(polynom)-1)
    else:
        p = proof([test[0], [[train[0][i], train[0][i+1]] for i in range(len(train[0])-1)] ] , polynom, len(polynom)-1)
    
    # Caculate the average absolute error
    m, s = error(test[1], p)
    print(jump, m, s, float(m+s), float(m-s))

    # Show the interpolation method graphic
    graphics(polynom, train, test, f"Pronostico del tiempo - Jump = {jump}")
    
main()