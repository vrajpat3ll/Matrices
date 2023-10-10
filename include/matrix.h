#ifndef MATRIX_H_
#define MATRIX_H_ 1

#include <bits/stdc++.h>
using namespace std;

class matrix
{
    int rows, columns;
    vector<vector<double>> M;

public:
    matrix();
    matrix(int m, int n);
    matrix(int n);

    /// @brief vector-like functionality
    /// @param i index
    /// @return i-th row of matrix
    vector<double> &operator[](int i) noexcept { return this->M[i]; }
    friend istream &operator>>(istream &in, matrix &a);
    friend ostream &operator<<(ostream &out, matrix &a);
    bool operator==(matrix &a);
    bool operator!=(matrix &a);
    matrix operator+(matrix &a);
    matrix operator-(matrix &a);
    matrix operator*(const matrix &a) const;
    matrix operator/(double divisor);

    matrix operator+=(matrix &a) { return (*this + a); }
    matrix operator-=(matrix &a) { return (*this - a); }
    matrix operator*=(matrix &a) { return (*this * a); }
    matrix operator^(int a);

    bool isSquare();
    bool isRowvector();
    bool isColumnvector();
    bool isDiagonal();
    friend pair<matrix, matrix> LUdecompose_crout(matrix &a);
    friend pair<matrix, matrix> LUdecompose_dolittle(matrix &a);
    friend matrix forward_sweep(matrix &, matrix &);
    friend matrix backward_sweep(matrix &, matrix &);
    friend matrix solveAXB(matrix &, matrix &);
    matrix identity(int size);
    // not complete
    // friend matrix Thomas_for_tri_diagonal(matrix &a, matrix &b);
    friend matrix solveAXB_gauss(matrix &A, matrix &B);
    friend matrix transpose(matrix &A);
    double determinant();
    matrix cofactormatrix(int row, int col);
};

/// @brief default constructor (3X3 matrix)
matrix::matrix()
{
    rows = 3;
    columns = 3;
    this->M = vector<vector<double>>(rows, vector<double>(columns));
}

/// @brief custom constructor (mXn matrix)
/// @param m number of rows
/// @param n number of columns
matrix::matrix(int m, int n)
{
    rows = m;
    columns = n;
    this->M = vector<vector<double>>(rows, vector<double>(columns));
}

/// @brief initialise a nXn square matrix
/// @param n
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

/// @brief decomposes a square matrix into a pair of lower-triangular and upper-triangular matrix using Crout's Decomposition
/// such that L.U = Original Matrix
/// @return L in first and U in second
pair<matrix, matrix> LUdecompose_crout(matrix &a)
{
    if (!a.isSquare())
    {
        cerr << "Cannot decompose." << endl;
        return {a, a};
    }

    matrix L(a.rows), U(a.columns);
    // Initialisation of L & U
    for (int i = 0; i < a.rows; i++)
    {
        for (int j = 0; j < a.columns; j++)
        {
            L[i][j] = 0;
            U[i][j] = (i == j) ? 1 : 0;
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

/// NOT FINISHED! <----<>---->
/// @brief decomposes a square matrix into a pair of lower-triangular and upper-triangular matrix using Doliitle Decomposition
/// such that L.U = Original Matrix
/// @return L in first and U in second
pair<matrix, matrix> LUdecompose_dolittle(matrix &a)
{
    if (!a.isSquare())
    {
        cerr << "Cannot decompose." << endl;
        return {a, a};
    }

    matrix L(a.rows), U(a.columns);

    // Initialisation of L & U
    for (int i = 0; i < a.rows; i++)
    {
        for (int j = 0; j < a.columns; j++)
        {
            L[i][j] = (i == j) ? 1 : 0;
            U[i][j] = 0;
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
    double x;
    cout << "Enter number of rows: ";
    in >> a.rows;
    cout << "Enter number of columns: ";
    in >> a.columns;
    a.M = vector<vector<double>>(a.rows, vector<double>(a.columns));

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
            out << (a[i][j] == -0.0) ? 0 : a[i][j];
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

// works
matrix forward_sweep(matrix &A, matrix &B)
{
    if (!B.isColumnvector() || A.rows != B.rows)
    {
        cout << "Cannot evaluate: Dimensions invalid." << endl;
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

// works
matrix backward_sweep(matrix &A, matrix &B)
{
    if (!B.isColumnvector() || A.rows != B.rows)
    {
        cout << "Cannot evaluate: Dimensions invalid." << endl;
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
    pair<matrix, matrix> LU = LUdecompose_crout(A);

    d = forward_sweep(LU.first, B);
    X = backward_sweep(LU.second, d);
    return X;
}

matrix solveAXB_gauss(matrix &A, matrix &B)
{
    if (!A.isSquare() || !B.isColumnvector() || A.rows != B.rows)
    {
        cerr << "Cannot solve the system: Invalid dimensions." << endl;
        return A; // Return the original matrix as an error indicator.
    }
    matrix a = A, b = B;
    int n = A.rows;
    matrix X(n, 1);

    // Forward elimination
    for (int i = 0; i < n - 1; i++)
    {
        for (int k = i + 1; k < n; k++)
        {
            double factor = A[k][i] / A[i][i]; // issue occurs when diagonal element is 0
            for (int j = i; j < n; j++)        // can be solved if we use pivot finding
            {
                a[k][j] -= factor * a[i][j];
            }
            b[k][0] -= factor * b[i][0];
        }
    }

    // Backward substitution
    for (int i = n - 1; i >= 0; i--)
    {
        double sum = 0;
        for (int j = i + 1; j < n; j++)
        {
            sum += a[i][j] * X[j][0];
        }
        X[i][0] = (b[i][0] - sum) / a[i][i];
    }

    return X;
}

bool matrix::operator==(matrix &a)
{
    if (a.rows != this->rows || a.columns != this->columns)
        return false;

    for (int i = 0; i < a.rows; i++)
    {
        for (int j = 0; j < a.columns; j++)
        {
            if (a[i][j] != this->M[i][j])
                return false;
        }
    }
    return true;
}

bool matrix::operator!=(matrix &b)
{
    if (*this == b)
        return false;
    return true;
}

matrix matrix::operator+(matrix &a)
{
    if (this->rows != a.rows || this->columns != a.columns)
    {
        cerr << "Cannot evaluate. Invalid dimensions." << endl;
        return *this;
    }
    matrix b(this->rows, this->columns);
    for (int i = 0; i < this->rows; i++)
    {
        for (int j = 0; j < this->columns; j++)
        {
            b[i][j] = this->M[i][j] + a[i][j];
        }
    }
    return b;
}

matrix matrix::operator-(matrix &a)
{
    if (this->rows != a.rows || this->columns != a.columns)
    {
        cerr << "Cannot evaluate. Invalid dimensions." << endl;
        return *this;
    }
    matrix b(this->rows, this->columns);
    for (int i = 0; i < this->rows; i++)
    {
        for (int j = 0; j < this->columns; j++)
        {
            b[i][j] = this->M[i][j] - a[i][j];
        }
    }
    return b;
}

matrix matrix::operator*(const matrix &a) const
{
    if (this->columns != a.rows)
    {
        cerr << "Matrix multiplication not possible: Invalid dimensions." << endl;
        return *this; // Return the original matrix as an error indicator.
    }

    matrix result(this->rows, a.columns);

    for (int i = 0; i < this->rows; i++)
    {
        for (int j = 0; j < a.columns; j++)
        {
            result[i][j] = 0;
            for (int k = 0; k < this->columns; k++)
            {
                result[i][j] += this->M[i][k] * a.M[k][j];
            }
        }
    }

    return result;
}

matrix matrix::operator/(double divisor)
{
    for (int i = 0; i < this->rows; i++)
    {
        for (int j = 0; j < this->columns; j++)
        {
            this->M[i][j] = this->M[i][j] / divisor;
        }
    }
    return *this;
}

matrix matrix::identity(int size)
{
    matrix b(size);
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            b[i][j] = (i == j) ? 1 : 0;
        }
    }
    return b;
}

matrix matrix::operator^(int a)
{
    if (!this->isSquare())
    {
        cerr << "Cannot evaluate. Invalid dimensions." << endl;
        return *this;
    }
    matrix b = identity(this->rows);

    for (int i = 0; i < a; i++)
    {
        b *= (*this);
    }
    return b;
}

double matrix::determinant()
{
    if (!this->isSquare())
    {
        cerr << "Determinant calculation is only valid for square matrices." << endl;
        return 0.0;
    }

    int n = this->rows;
    switch (n)
    {
    case 1:
        return this->M[0][0];
    case 2:
        return this->M[0][0] * this->M[1][1] - this->M[1][0] * this->M[0][1];

    default:
        break;
    }
    double det = 0;
    for (int i = 0; i < n; i++)
    {
        double cofactor = this->M[0][i] * this->cofactormatrix(0, i).determinant();
        det += (i & 1) ? cofactor : -cofactor;
    }

    return det;
}

matrix matrix::cofactormatrix(int row, int col)

{
    int n = this->rows;
    matrix cofactor(n - 1, n - 1);
    int newRow = 0, newCol = 0;

    for (int i = 0; i < n; i++)
    {
        if (i == row)
            continue;

        newCol = 0;
        for (int j = 0; j < n; j++)
        {
            if (j == col)
                continue;

            cofactor[newRow][newCol] = this->M[i][j];
            newCol++;
        }
        newRow++;
    }

    return cofactor;
}

matrix transpose(matrix &A)
{
    matrix B(A.columns, A.rows);
    for (int i = 0; i < A.rows; i++)
    {
        for (int j = 0; j < A.columns; j++)
        {
            B[j][i] = A[i][j];
        }
    }
    return B;
}
#endif