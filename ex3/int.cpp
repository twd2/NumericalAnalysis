#include <cmath>
#include <iostream>
#include <iomanip>
#include <functional>

// 梯形公式
double T(std::function<double (double)> &f, double a, double b, double n)
{
    double h = (b - a) / n;
    double middle = 0;
    double x = a + h;
    for (int i = 1; i <= n - 1; ++i)
    {
        middle += f(x);
        x += h;
    }
    return (f(a) + f(b) + 2 * middle) * h / 2;
}

// 辛普森公式
double S(std::function<double (double)> &f, double a, double b, double n)
{
    double h = (b - a) / n;
    double middle1 = 0;
    double x = a + 0.5 * h;
    for (int i = 0; i <= n - 1; ++i)
    {
        middle1 += f(x);
        x += h;
    }
    double middle2 = 0;
    x = a + h;
    for (int i = 1; i <= n - 1; ++i)
    {
        middle2 += f(x);
        x += h;
    }
    return (f(a) + f(b) + 2 * middle2 + 4 * middle1) * h / 6;
}

// 龙贝格算法
double L(std::function<double (double)> &f, double a, double b, double n)
{
    // TODO
    return 0.0;
}

int main()
{
    std::cout << std::fixed << std::setprecision(20);
    std::function<double (double)> f = [](double x) { return exp(x); };
    double acc_result = exp(1) - 1;
    std::cout << "e - 1 = " << acc_result << std::endl;
    double t_result = T(f, 0.0, 1.0, 1000);
    std::cout << "  T   = " << t_result << ", diff = " << t_result - acc_result << std::endl;
    double s_result = S(f, 0.0, 1.0, 8);
    std::cout << "  S   = " << s_result << ", diff = " << s_result - acc_result << std::endl;
    double l_result = L(f, 0.0, 1.0, 5);
    std::cout << "  L   = " << l_result << ", diff = " << l_result - acc_result << std::endl;

    int n = 1;
    double eps = 1e-6;
    double result;
    do
    {
        result = T(f, 0.0, 1.0, n);
        std::cout << "n = " << n << ", " << "T = " << result << ", diff = " << result - acc_result << std::endl;
        ++n;
    }
    while (fabs(result - acc_result) >= eps);
    n = 1;
    do
    {
        result = S(f, 0.0, 1.0, n);
        std::cout << "n = " << n << ", " << "S = " << result << ", diff = " << result - acc_result << std::endl;
        ++n;
    }
    while (fabs(result - acc_result) >= eps);
    /*n = 1;
    do
    {
        result = L(f, 0.0, 1.0, n);
        std::cout << "n = " << n << ", " << "L = " << result << ", diff = " << result - acc_result << std::endl;
        ++n;
    }
    while (fabs(result - acc_result) >= eps);*/
    return 0;
}