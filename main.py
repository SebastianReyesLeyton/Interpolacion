
from numpy import mean, std
from pandas import read_csv, DataFrame
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from interpolations import (
    polynomial
)
from functions import (
    graphics
)
from time import time

def main():

    # Upload data of .csv file
    data = read_csv("./DataSet1.csv", header=None) # The time is expressed on days.

    # Obtain the array of time and its linked values
    T = data.drop(1, axis=1)
    y = data.drop(0, axis=1)

    # Obtain the first value of time
    d = T[0][0]

    # Transform the time array
    for i in range(T.shape[0]): T[0][i] -= d

    print(T[0])

    # Create the subset of train and test, and their corresponds outputs
    T_train, T_test, y_train, y_test = train_test_split(T, y, test_size=0.97, random_state=42)

    # Create the train matrix
    train = [list(T_train[0]), list(y_train[1])]

    print(len(T_train))

    # Create the test matrix
    test = [list(T_test[0]), list(y_test[1])]

    Data = [[-1, -0.5, 0.0, 0.5, 1], [1, 0.5, 0.0, 0.5, 2.0]]

    # Obtain the polynomial through normalize equations
    start = time()
    polynom = polynomial(test)
    stop = time()

    print(polynom)

    print(f"Time: {stop - start}")

    graphics(polynom, train, "Ejemplo")
    graphics(polynom, test, "Ejemplo")

main()