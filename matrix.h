#ifndef MATRIX_H_
#define MATRIX_H_ 1

#include <bits/stdc++.h>

using namespace std;

/// @brief wraps rows,columns & 2D array
class matrix
{
    int rows, columns;
    vector<vector<double>> M;

public:
    /// @brief default constructor (3X3 matrix)
    matrix();

    /// @brief custom constructor (mXn matrix)
    /// @param m number of rows
    /// @param n number of columns
    matrix(int m, int n);

    /// @brief initialise a nXn square matrix
    /// @param n
    matrix(int n);

    /// @brief vector-like functionality
    /// @param i index
    /// @return i-th row of matrix
    vector<double> &operator[](int i) noexcept { return this->M[i]; }

    friend istream &operator>>(istream &in, matrix &a);
    friend ostream &operator<<(ostream &out, matrix &a);

    friend pair<matrix, matrix> LUdecompose(matrix &a);

    bool isSquare();
    bool isRowvector();
    bool isColumnvector();
    bool isDiagonal();

    friend matrix forward_sweep(matrix &, matrix &);
    friend matrix backward_sweep(matrix &, matrix &);
    friend matrix solveAXB(matrix &, matrix &);

    friend bool operator==(matrix &a, matrix &b);
    friend bool operator!=(matrix &a, matrix &b);
};

/// @brief decomposes a square matrix into a pair of lower-triangular(L) and upper-triangular matrix(U)
/// @return L in first and U in second
pair<matrix, matrix> LUdecompose(matrix &a)
{
    if (!a.isSquare())
    {
        cerr << "Cannot decompose." << endl;
        return {a, a}; // I didn;t know what to return, so i just return matrix a to it
    }

    matrix L(a.rows), U(a.columns);
    // Initialisation of L & U
    for (int i = 0; i < a.rows; i++)
    {
        for (int j = 0; j < a.columns; j++)
        {
            L[i][j] = 0;
            if (i == j)
            {
                U[i][j] = 1;
            }
            else
            {
                U[i][j] = 0;
            }
        }
    }

    // Decompostion Operations
    for (int i = 0; i < a.rows; i++)
    {
        L[i][0] = a[i][0];
    }
    for (int i = 1; i < a.columns; i++)
    {
        U[0][i] = a[0][i] / L[0][0];
    }

    // L operation
    for (int j = 2; j <= a.rows - 1; j++)
    {
        for (int i = j; i <= a.columns; i++)
        {
            L[i - 1][j - 1] = a[i - 1][j - 1];
            for (int k = 1; k <= j - 1; k++)
            {
                L[i - 1][j - 1] = L[i - 1][j - 1] - L[i - 1][k - 1] * U[k - 1][j - 1];
            }
        }
        // U operation
        for (int x = j + 1; x <= a.rows; x++)
        {
            U[j - 1][x - 1] = a[j - 1][x - 1];
            for (int w = 1; w <= j - 1; w++)
            {
                U[j - 1][x - 1] -= L[j - 1][w - 1] * U[w - 1][x - 1];
            }
            U[j - 1][x - 1] /= L[j - 1][j - 1];
        }
    }

    // L[n-1][n-1] operation
    L[a.rows - 1][a.columns - 1] = a[a.rows - 1][a.columns - 1];
    for (int i = 1; i <= a.rows - 1; i++)
    {
        L[a.rows - 1][a.columns - 1] -= L[a.rows - 1][i - 1] * U[i - 1][a.columns - 1];
    }

    return {L, U};
}

/// @brief takes input of a matrix
/// @param a matrix
istream &operator>>(istream &in, matrix &a)
{
    do
    {
        cout << "Enter number of rows: ";
        in >> a.rows;
    } while (a.rows <= 0);

    do
    {
        cout << "Enter number of columns: ";
        in >> a.columns;
    } while (a.columns <= 0);

    a.M = vector<vector<double>>(a.rows, vector<double>(a.columns));

    double x = 0;
    for (int i = 0; i < a.rows; i++)
    {
        for (int j = 0; j < a.columns; j++)
        {
            cout << "Enter value at (" << i + 1 << "," << j + 1 << ") : ";
            in >> x;
            a.M[i][j] = x;
        }
    }

    return in;
}

/// @brief prints the matrix and also the rows and columns
/// @param a matrix
ostream &operator<<(ostream &out, matrix &a)
{
    out << "Rows: " << a.rows << endl;
    out << "Columns: " << a.columns << endl;
    out << "[";
    for (int i = 0; i < a.rows; i++)
    {
        out << "[";
        for (int j = 0; j < a.columns; j++)
        {
            out << a.M[i][j];
            if (j != a.columns - 1)
                out << ",";
        }
        out << "]";
        if (i != a.rows - 1)
        {
            out << ",";
            out << endl;
        }
    }
    out << "]" << endl;

    return out;
}

matrix forward_sweep(matrix &A, matrix &B)
{
    if (!B.isColumnvector() || A.rows != B.rows)
    {
        cout << "Cannot evaluate." << endl;
        return A;
    }
    matrix X(A.rows, 1);
    X[0][0] = B[0][0] / A[0][0];

    for (int i = 2; i <= A.rows; i++)
    {
        int sum = 0;
        for (int j = 1; j <= i - 1; j++)
        {
            sum = sum + A[i - 1][j - 1] * X[j - 1][0];
        }
        X[i - 1][0] = (B[i - 1][0] - sum) / A[i - 1][i - 1];
    }
    return X;
}

matrix backward_sweep(matrix &A, matrix &B)
{
    if (!B.isColumnvector() || A.rows != B.rows)
    {
        cout << "Cannot evaluate." << endl;
        return A;
    }

    matrix X(A.rows, 1);
    X[A.rows - 1][0] = B[A.rows - 1][0] / A[A.rows - 1][A.rows - 1];

    for (int i = A.rows - 1; i >= 1; i--)
    {
        int sum = 0;
        for (int j = i + 1; j <= A.rows; j++)
        {
            sum = sum + A[i - 1][j - 1] * X[j - 1][0];
        }
        X[i - 1][0] = (B[i - 1][0] - sum) / A[i - 1][i - 1];
    }
    return X;
}

/// @brief solves the equation A.X = B using LU decompostion
/// @param A Coefficient Matrix
/// @param B Source Vector
/// @return Solution Set of the equation
matrix solveAXB(matrix &A, matrix &B)
{
    if (!B.isColumnvector() || A.rows != B.rows)
    {
        cerr << "Cannot compute." << endl;
        return A;
    }
    matrix d, X;
    pair<matrix, matrix> LU = LUdecompose(A);

    d = forward_sweep(LU.first, B);
    X = backward_sweep(LU.second, d);
    return X;
}

bool operator==(matrix &a, matrix &b)
{
    if (a.rows != b.rows || a.columns != b.columns)
        return false;

    for (int i = 0; i < a.rows; i++)
    {
        for (int j = 0; j < a.columns; j++)
        {
            if (a[i][j] != b[i][j])
                return false;
        }
    }
    return true;
}

bool operator!=(matrix &a, matrix &b)
{
    if (a == b)
        return false;
    return true;
}

matrix::matrix()
{
    rows = 3;
    columns = 3;
    this->M = vector<vector<double>>(rows, vector<double>(columns));
}

matrix::matrix(int m, int n)
{
    rows = m;
    columns = n;
    this->M = vector<vector<double>>(rows, vector<double>(columns));
}

matrix::matrix(int n)
{
    rows = n;
    columns = n;
    this->M = vector<vector<double>>(rows, vector<double>(columns));
}

bool matrix::isSquare()
{
    if (this->rows != this->columns)
        return false;
    return true;
}

bool matrix::isRowvector()
{
    if (this->columns == 1)
        return true;
    return false;
}

bool matrix::isColumnvector()
{
    if (this->columns == 1)
        return true;
    return false;
}

bool matrix::isDiagonal()
{
    for (int i = 0; i < this->rows; i++)
    {
        for (int j = 0; j < this->columns; j++)
        {
            if (i != j && this->M[i][j] != 0)
            {
                return false;
            }
        }
    }
    return true;
}
#endif