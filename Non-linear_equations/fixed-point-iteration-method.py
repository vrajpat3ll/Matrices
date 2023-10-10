# for use of mathematical constants
import numpy as np

TOLERANCE = 1e-6
MAX_ITER = 1000


def f(x):  # this would be given
    return x**2 - 2 * x + 1


def g(x: float):  # have to be found out
    return float((2 * x - 1) ** 0.5)


def fixd_pt_iter_method():
    x = float(input("Enter guess for root of function: "))

    print("-"*129)
    print("Iteration\t|\tx\t\t|\tg(x)\t\t|\tf(x)\t\t\t|")

    for i in range(1, MAX_ITER + 1):
        print("-"*129)
        print(i, "\t\t|\t", format(x, ".6f"), "\t|\t", format(
            g(x), ".6f"), "\t|\t", format(f(x), ".7f"), "\t\t|", sep="")

        x = g(x)  # algorithm

        if (abs(f(x)) <= TOLERANCE):
            print("-"*129)
            print(i + 1, "\t\t|\t", format(x, ".6f"), "\t|\t", format(g(x),
                  ".6f"), "\t|\t", format(f(x)), "\t|", sep="")
            print("-"*129)
            break
    return


if __name__ == '__main__':
    fixd_pt_iter_method()
