// TESTING NOT COMPLETED!

#include "../../include/matrix.h"

#define MAX_ITER 50

int main()
{
    int N_Unknowns, N_Equations;
    cout << "Enter the number of unknowns in the system: ";
    cin >> N_Unknowns;
    cout << "Enter the number of equations in the system: ";
    cin >> N_Equations;
    if (N_Equations > N_Unknowns)
    {
        cout << "Can't solve the system. More equations than unknowns." << endl;
        return 1;
    }
    if (N_Equations < N_Unknowns)
    {
        cout << "Can't solve the system. More unknowns than equations." << endl;
        return -1;
    }

    matrix Coefficients(N_Equations, N_Unknowns), source(N_Unknowns, 1), X(N_Unknowns, 1);

    cout << "\t\tSolve AX=B" << endl
         << endl;
    cout << "Enter the coffients-matrix: A" << endl;
    cin >> Coefficients;
    cout << "Enter source-vector: B" << endl;
    cin >> source;
    cout << "Enter initial guess-vector: X" << endl;
    cin >> X;
    matrix x = X;

    // Jacobi Algorithm
    for (int iter = 0; iter < MAX_ITER; iter++)
    {
        for (int i = 0; i < N_Equations; i++)
        {
            X[i][0] = x[i][0] + (source[i][0] - (Coefficients * x)[i][0]) / Coefficients[i][i];
        }

        X = x;
    }
    cout << "Solution: X is" << endl;
    cout << X << endl;
}