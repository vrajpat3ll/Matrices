#ifndef INTERPOLATION_H_
#define INTERPOLATION_H_
#include "matrix.h"

std::pair<double, std::string> S(int i, double x)
{
    switch (i)
    {
    case 0:
        return {1, ""};
    case 1:
        return {x, "x"};
    case 2:
        return {pow(x, 2), "x^2"};
    case 3:
        return {pow(x, 3), "x^3"};
    case 4:
        return {pow(x, 4), "x^4"};
    case 5:
        return {pow(x, 5), "x^5"};
    case 6:
        return {pow(x, 6), "x^6"};
    case 7:
        return {pow(x, 7), "x^7"};
    case 8:
        return {pow(x, 8), "x^8"};
    case 9:
        return {pow(x, 9), "x^9"};
    case 10:
        return {pow(x, 10), "x^10"};
    default:
        break;
    }
    return {0, ""};
}

matrix linearRegressionCoefficients(std::pair<double, std::string> (*func)(int, double), std::vector<std::pair<double, double>> points, int N_degree)
{
    int N_points = points.size();
    matrix F(N_points, N_degree);
    for (int i = 0; i < N_points; i++)
        for (int j = 0; j < N_degree; j++)
            F[i][j] = func(j, points[i].first).first;

    matrix Ft = transpose(F);

    // now we do A=Ft*F which is equal to points.second
    matrix A = Ft * F;

    matrix B(N_points, 1);
    for (int i = 0; i < N_points; i++)
        B[i][0] = points[i].second; // y-values

    matrix C = Ft * B;
    matrix regressionCoefficients = solveAXB_gauss(A, C);

    return regressionCoefficients;
}

void printRegressionFunction(matrix regressionCoefficients, std::pair<double, std::string> (*func)(int, double), int N_degree)
{
    std::cout << "f(x) = ";
    for (int i = 0; i < N_degree; i++)
    {                                                                   // func(i,0) me 0 se koi farak nahi padta
        std::cout << regressionCoefficients[i][0] << func(i, 0).second; // used func(i,0) here. bcoz we only want std::string part of func()
        if (i != N_degree - 1)                                          // function and 0 does not play any role here.
            std::cout << " + ";
    }
    std::cout << std::endl;
}

double interpolateRegression(std::pair<double, std::string> (*func)(int, double), matrix &regressionCoefficients, double x_val)
{
    double f_x = 0;
    for (int i = 0; i < regressionCoefficients.get_rows(); i++)
        f_x += regressionCoefficients[i][0] * func(i, x_val).first;
    return f_x;
}

// helper for nth-order lagrange interpolation
double L_i(int i, double x, int N_order, std::vector<double> points)
{
    double l = 1;
    for (int k = 0; k < N_order + 1; k++)
    {
        if (k == i)
            continue;
        l *= (x - points[k]) / (points[i] - points[k]);
    }
    return l;
}

double firstOrderLagrangeInterpolation(std::vector<double> f_x, std::vector<double> points, double x)
{
    double value = 0;
    value += (x - points[1]) / (points[0] - points[1]) * f_x[0];
    value += (x - points[0]) / (points[1] - points[0]) * f_x[1];

    return value;
}
double secondOrderLagrangeInterpolation(std::vector<double> f_x, std::vector<double> points, double x)
{
    double value = 0;
    value += (x - points[1]) * (x - points[2]) * f_x[0] / ((points[0] - points[1]) * (points[0] - points[2]));
    value += (x - points[0]) * (x - points[2]) * f_x[1] / ((points[1] - points[0]) * (points[1] - points[2]));
    value += (x - points[0]) * (x - points[1]) * f_x[2] / ((points[2] - points[0]) * (points[2] - points[1]));
    return value;
}


/// @param f_x value at point[i] 
/// @param points point[i]
/// @param x point at which value is to be evaluated
/// @param N_order order of lagrange interpolation
/// @return interpolated value
double lagrangeInterpolation(std::vector<double> f_x, std::vector<double> points, double x, int N_order)
{
    double value = 0;
    for (int i = 0; i <= N_order; i++)
        value += L_i(i, x, N_order, points) * f_x[i];
    return value;
}
//TODO: Spline interpolation implementation
double splineInterpolation{};

#endif