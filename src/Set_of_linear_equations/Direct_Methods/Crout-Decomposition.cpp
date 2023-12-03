#include "../../include/matrix.h"

int main()
{
    matrix A;
    cin >> A;
    pair<matrix,matrix> LU = LUdecompose_crout(A);
    matrix Lower = LU.first;
    matrix Upper = LU.second;
    cout << Lower << endl;
    cout << Upper << endl;
}