#include <cmath>
#include <iostream>
#include <iomanip>
#include <functional>

// Gauss公式
double G(std::function<double (double)> &f, double a, double b, double n)
{
    double h = (b - a) / n;
    double sum = 0;
    double x1 = a + 0.5 * h - h * sqrt(3) / 6, x2 = a + 0.5 * h + h * sqrt(3) / 6;
    for (int i = 0; i <= n - 1; ++i)
    {
        sum += f(x1) + f(x2);
        x1 += h;
        x2 += h;
    }
    return sum * h / 2;
}

int main()
{
    std::cout << std::fixed << std::setprecision(20);
    std::function<double (double)> f = [](double x) { return 4 / (1 + x * x); };
    double acc_result = M_PI;
    std::cout << "pi = " << acc_result << std::endl;

    double g_result;
    g_result = G(f, 0.0, 1.0, 13);
    std::cout << "G = " << g_result << ", diff = " << g_result - acc_result << std::endl;
    return 0;
}