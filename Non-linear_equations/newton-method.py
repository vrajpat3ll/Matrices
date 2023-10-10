import numpy as np

TOLERANCE = 1e-6
MAX_ITER = 1000
DELX = 1e-4


def f(x):
    return x**2 - 2 * x + 1


def derivative(f, x0: float):
    # using the definition of derivative of f(x) at x0
    return (f(x0 + DELX) - f(x0)) / DELX


def newtons_method():
    init_guess = float(input("Enter guess for root of function: "))

    print("-"*113)
    print("Iteration\t|\tRoot\t\t|\tf(x)\t\t\t|\tf'(x)\t\t\t\t|")

    for i in range(1, MAX_ITER + 1):

        Xnew = init_guess - f(init_guess) / derivative(f, init_guess)
        yn = f(Xnew)

        print("-"*113)
        print(i, "\t\t|\t", format(init_guess, ".6f"),
              "\t|\t", format(f(init_guess), ".7f"), "\t\t|\t", (derivative(f, init_guess)), "\t\t|", sep="")

        init_guess = Xnew

        if (abs(yn) <= TOLERANCE):
            print("-"*113)
            print(i + 1, "\t\t|\t", format(init_guess, ".6f"),
                  "\t|\t", format(f(init_guess)), "\t|\t", derivative(f, init_guess), "\t\t|", sep="")
            print("-"*113)

            break
    return


if __name__ == '__main__':
    newtons_method()
