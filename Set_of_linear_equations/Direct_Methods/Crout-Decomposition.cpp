#include "../../include/matrix.h"

int main()
{
    matrix A;
    cin >> A;
    matrix Lower = LUdecompose_crout(A).first;
    matrix Upper = LUdecompose_crout(A).second;
    cout << Lower << endl;
    cout << Upper << endl;
}