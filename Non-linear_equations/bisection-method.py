MAX_ITER = 1000
DELX = 1
TOLERANCE = 1e-12


def bracket_roots(f):
    x = -1e2
    while f(x) * f(x + DELX) > 0:
        x = x + DELX
    return x, x + DELX


def f(x):
    return x**2 - 2 * x + 1


def bisection_method():
    xp, xn = bracket_roots(f)
    print("-" * 81)
    print("Iteration\t|\tRoot\t\t|\tf(x)\t\t\t\t|")
    for i in range(1, MAX_ITER + 1):
        Xnew = (xp + xn) / 2
        y = f(Xnew)

        print("-" * 81)
        print(i, "\t\t|\t", format(Xnew, ".6f"),"\t|\t", format(y), "\t\t|", sep="")
        
        if abs(y) <= TOLERANCE:
            print("-" * 81)
            return
        
        if f(xp) < 0:
            if y < 0:
                xp = Xnew
            else:
                xn = Xnew

        else:
            if y < 0:
                xn = Xnew
            else:
                xp = Xnew

    return


if __name__ == '__main__':
    bisection_method()
