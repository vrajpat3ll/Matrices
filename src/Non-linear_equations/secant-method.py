# for use of mathematical constants
import numpy as np

TOLERANCE = 1e-6
MAX_ITER = 1000


def f(x):
    return x**2 - 2 * x + 1


def secant_method():
    x1 = float(input("Enter first guess for root of function: "))
    x2 = float(input("Enter second guess for root of function: "))

    print("-"*129)
    print("Iteration\t|\tx1\t\t|\tx2\t\t|\tf(x2)\t\t\t|\tslope\t\t\t|")

    for i in range(1, MAX_ITER + 1):
        m = (f(x2) - f(x1)) / (x2 - x1)
        Xnew = x2 - f(x2) / m
        yn = f(Xnew)

        print("-"*129)
        print(i, "\t\t|\t", format(x1, ".6f"), "\t|\t", format(x2, ".6f"),
              "\t|\t", format(f(x2), ".7f"), "\t\t|\t", m, "\t|", sep="")
        x1 = x2
        x2 = Xnew

        if (abs(yn) <= TOLERANCE):
            print("-"*129)
            print(i + 1, "\t\t|\t", format(x1, ".6f"), "\t|\t", format(x2,
                  ".6f"), "\t|\t", format(f(x2)), "\t|\t", m, "\t|", sep="")
            print("-"*129)
            break
    return


if __name__ == '__main__':
    secant_method()