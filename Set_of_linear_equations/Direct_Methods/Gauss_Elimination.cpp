#include "../../include/matrix.h"

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

    matrix Coefficients(N_Equations, N_Unknowns), source(N_Unknowns, 1);
    cout << "\t\tSolve AX=B" << endl;
    cout << "Enter the coffients-matrix: A" << endl;
    cin >> Coefficients;
    cout << "Enter source-vector: B" << endl;
    cin >> source;
    matrix X = solveAXB_gauss(Coefficients, source);
    cout << "" << endl;
    cout << X;
}