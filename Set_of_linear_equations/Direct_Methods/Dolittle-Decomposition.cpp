#include "../../include/matrix.h"

int main()
{
    matrix A;
    cin >> A;
    matrix Lower = LUdecompose_dolittle(A).first;
    matrix Upper = LUdecompose_dolittle(A).second;
    cout << Lower << endl;
    cout << Upper << endl;
}