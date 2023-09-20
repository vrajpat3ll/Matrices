#include "../matrix.h"

const double h = 2e-3;
const double x10 = 1;
const double k1 = 1, k2 = 0.2, k3 = 0.05, k4 = 0.4;
const int MAX_ITER = 100;

double f(int i, double x1, double x2, double x3, double x4, double x5);

double partial_derivative(double (*func)(int, double, double, double, double, double), int i, int j, double x1, double x2, double x3, double x4, double x5);

vector<double> &solveEQns(double (*f)(int, double, double, double, double, double), vector<double> &guess, int N_unknowns);

template <class T>
ostream &operator<<(ostream &out, vector<T> &v);

// MAIN FUNCTION
int main()
{
    vector<double> guess(4, 1);
    cout << "Iteration\tRoot\t\t|\tDX" << endl;
    solveEQns(f, guess, 4);
    cout << guess;
}

/// @brief system of (usually) non-linear equations ( f(x_i,...) = 0 )
/// @param i index of expression (from 0)
/// @param x1 variable 1
/// @param x2 variable 2
/// @param x3 variable 3
/// @param x4 variable 4
/// @param x5 variable 5
/// @return value of expression
double f(int i, double x1, double x2, double x3, double x4, double x5)
{
    switch (i)
    {
    case 0:
        return -x1 + x10 + 2 * (-k1 * x1 - k2 * pow(x1, 1.5) + k3 * pow(x3, 2));
    case 1:
        return -x2 + 2 * (2 * k1 * x1 - k4 * pow(x2, 2));
    case 2:
        return -x3 + 2 * (k2 * pow(x1, 1.5) + k4 * pow(x2, 2) - k3 * pow(x3, 2));
    case 3:
        return -x4 + 2 * (k4 * pow(x2, 2));
    default:
        return 0;
    }
    return 0;
}

/// @brief partial derivative approximation ( ∆x = 2e-3 )
/// @param func function pointer of a function
/// @param i index of expression
/// @param j index of variable
/// @param x1 variable 1
/// @param x2 variable 2
/// @param x3 variable 3
/// @param x4 variable 4
/// @param x5 variable 5
/// @return df_i/dx_j
double partial_derivative(double (*func)(int, double, double, double, double, double), int i, int j, double x1, double x2, double x3, double x4, double x5)
{
    switch (j)
    {
    case 0:
        return (func(i, x1 * (1 + h), x2, x3, x4, x5) - func(i, x1, x2, x3, x4, x5)) / (h * x1);
    case 1:
        return (func(i, x1, x2 * (1 + h), x3, x4, x5) - func(i, x1, x2, x3, x4, x5)) / (h * x2);
    case 2:
        return (func(i, x1, x2, x3 * (1 + h), x4, x5) - func(i, x1, x2, x3, x4, x5)) / (h * x3);
    case 3:
        return (func(i, x1, x2, x3, x4 * (1 + h), x5) - func(i, x1, x2, x3, x4, x5)) / (h * x4);
    case 4:
        return (func(i, x1, x2, x3, x4, x5 * (1 + h)) - func(i, x1, x2, x3, x4, x5)) / (h * x5);

    default:
        return 0;
    }
}

/// @brief Newton-Raphson's method
/// @param f expression array
/// @param guess initial guess for the solution
/// @param N_unknowns number of unknowns (upto 5)
/// @return root of the equations
vector<double> &solveEQns(double (*f)(int, double, double, double, double, double), vector<double> &guess, int N_unknowns)
{
    for (int iter = 1; iter <= MAX_ITER; iter++)
    {

        matrix Jacobian(N_unknowns, N_unknowns), source(N_unknowns, 1), DX(N_unknowns, 1);

        for (int i = 0; i < N_unknowns; i++)
        {
            for (int j = 0; j < N_unknowns; j++)
            {
                Jacobian[i][j] = partial_derivative(f, i, j, guess[0], guess[1], guess[2], guess[3], guess[4]);
            }
        }

        for (int i = 0; i < N_unknowns; i++)
        {
            source[i][0] = -1 * f(i, guess[0], guess[1], guess[2], guess[3], guess[4]);
        }

        DX = solveAXB(Jacobian, source);

        for (int i = 0; i < N_unknowns; i++)
        {
            guess[i] += DX[i][0];
        }

        cout << iter;

        for (int i = 0; i < N_unknowns; i++)
        {
            cout << "\t\tx_" << i + 1 << " = " << guess[i] << "\t|\t" << DX[i][0] << endl;
        }
        cout << "____________________________________________________________________________" << endl;
        cout << endl;
    }
    return guess;
}

/// @brief outputs the
/// @tparam T can be int,long,long long,float,double,long double
/// @param v vector
template <class T>
ostream &operator<<(ostream &out, vector<T> &v)
{
    for (int i = 0; i < v.size(); i++)
    {
        out << "x_" << i + 1 << " = " << setprecision(10) << v[i] << endl;
    }
    return out;
}