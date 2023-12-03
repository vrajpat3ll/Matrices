#define newton_raphson_constraint_
// #ifndef newton_raphson_
#include "../../include/matrix.h"
// #endif
using namespace std;

int MAX_ITERS = 1000;
const double h = 2e-3;
const double x10 = 1;
const double k1 = 1, k2 = 0.2, k3 = 0.05, k4 = 0.4;
const int MAX_ITER = 100;

double f(int i, double x1, double x2, double x3, double x4, double x5);

double partial_derivative(double (*func)(int, double, double, double, double, double), double deltaVar, int indexOfFunction, int indexOfVariable, double x1, double x2, double x3, double x4, double x5);

vector<double> &solveEQ(double (*f)(int, double, double, double, double, double), vector<double> &guess, vector<double> &root, int N_unknowns);

template <class T>
ostream &operator<<(ostream &out, vector<T> &v);

// MAIN FUNCTION
// int main()
// {
//     vector<double> guess(4, 1);
//     vector<double> root = {0.318866, 0.783884, 0.534982, 0.491579};
//     cout << "---------------------------------------------------------------------------------------------------" << endl;
//     cout << "Iteration\t|\tRoot\t\t\t|\tResidue\t\t|\tDX" << endl;
//     cout << "---------------------------------------------------------------------------------------------------" << endl;
//     solveEQ(f, guess, root, 4);
//     cout << guess;
// }

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
    return 0; // 0 is returned because no value is to be calculated at that index
}

/// @brief partial derivative approximation
/// @param func function pointer of a function
/// @param indexOfFunction index of expression
/// @param indexOfVariable index of variable
/// @param deltaVar change in variable for approximation (∆x)
/// @param x1 variable 1
/// @param x2 variable 2
/// @param x3 variable 3
/// @param x4 variable 4
/// @param x5 variable 5
/// @returns df_i/dx_j
double partial_derivative(double (*func)(int, double, double, double, double, double), double deltaVar, int indexOfFunction, int indexOfVariable, double x1, double x2, double x3, double x4, double x5)
{
    switch (indexOfVariable)
    {
    case 0:
        return (func(indexOfFunction, x1 * (1 + deltaVar), x2, x3, x4, x5) - func(indexOfFunction, x1, x2, x3, x4, x5)) / (deltaVar * x1);
    case 1:
        return (func(indexOfFunction, x1, x2 * (1 + deltaVar), x3, x4, x5) - func(indexOfFunction, x1, x2, x3, x4, x5)) / (deltaVar * x2);
    case 2:
        return (func(indexOfFunction, x1, x2, x3 * (1 + deltaVar), x4, x5) - func(indexOfFunction, x1, x2, x3, x4, x5)) / (deltaVar * x3);
    case 3:
        return (func(indexOfFunction, x1, x2, x3, x4 * (1 + deltaVar), x5) - func(indexOfFunction, x1, x2, x3, x4, x5)) / (deltaVar * x4);
    case 4:
        return (func(indexOfFunction, x1, x2, x3, x4, x5 * (1 + deltaVar)) - func(indexOfFunction, x1, x2, x3, x4, x5)) / (deltaVar * x5);

    default:
        return 0;
    }
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

/// @brief Newton-Raphson's method for solving  upto 5 non-linear equations
/// @param f expression array
/// @param guess initial guess-values for the solution
/// @param N_unknowns number of unknowns (upto 5)
/// @return root of the equations
vector<double> &solveEQ(double (*f)(int, double, double, double, double, double), vector<double> &guess, vector<double> &root, int N_unknowns)
{
    vector<double> residue(N_unknowns);
    for (int iter = 1; iter <= MAX_ITERS; iter++)
    {
        matrix Jacobian(N_unknowns, N_unknowns), source(N_unknowns, 1), DX(N_unknowns, 1);
        for (int i = 0; i < N_unknowns; i++)
        {
            for (int j = 0; j < N_unknowns; j++)
            {
                Jacobian[i][j] = partial_derivative(f, h, i, j, guess[0], guess[1], guess[2], guess[3], guess[4]);
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

        for (int i = 0; i < N_unknowns; i++)
        {
            residue[i] = guess[i] - root[i];
        }

        cout << iter;

        for (int i = 0; i < N_unknowns; i++)
        {
            cout << "\t\t|\tx_" << i + 1 << " = " << guess[i] << "\t\t|\t" << residue[i] << "\t|\t" << DX[i][0] << endl;
        }

        cout << "___________________________________________________________________________________________________" << endl;
        cout << endl;

        if (abs(DX[0][0]) <= 1e-4 && abs(DX[1][0]) <= 1e-4 && abs(DX[2][0]) <= 1e-4 && abs(DX[3][0]) <= 1e-4)
        {
            return guess;
        }
    }
    return guess;
}