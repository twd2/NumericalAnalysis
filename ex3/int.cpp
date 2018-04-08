#include <cmath>
#include <iostream>
#include <iomanip>
#include <functional>
#include <vector>
#include <map>

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
double L(std::function<double (double)> &f, double a, double b, double eps)
{
    std::cout << std::fixed << std::setprecision(8);
    std::map<int, double [2]> Tm;
    Tm[0][1] = (f(a) + f(b)) * (b - a) * 0.5;
    std::cout << "0 " << Tm[0][1] << std::endl;
    int k = 0, n = 1;
    double h = b - a;
    while (true)
    {
        Tm[0][0] = Tm[0][1];
        // 递推公式计算 T_0^{(k)}
        double diff = 0.0, x = a + 0.5 * h;
        for (int i = 0; i < n; ++i)
        {
            diff += f(x);
            x += h;
        }
        Tm[0][1] = Tm[0][0] * 0.5 + diff * h * 0.5;
        ++k;
        n *= 2;
        h /= 2.0;
        std::cout << k << " " << Tm[0][1] << " ";

        // 计算 T_m^{(k)}
        double coeff = 4.0;
        for (int m = 1; m <= k; ++m)
        {
            Tm[m][0] = Tm[m][1];
            Tm[m][1] = (coeff * Tm[m - 1][1] - Tm[m - 1][0]) / (coeff - 1.0);
            std::cout << Tm[m][1] << " ";
            coeff *= 4.0;
        }
        std::cout << std::endl;
        if (std::abs(Tm[k][1] - Tm[k - 1][1]) < eps)
        {
            break;
        }
    }
    std::cout << std::fixed << std::setprecision(20);
    return Tm[k][1];
}

int main()
{
    const double eps = 1e-6;
    std::cout << std::fixed << std::setprecision(20);
    std::function<double (double)> f = [](double x) { return exp(x); };
    double acc_result = exp(1) - 1;
    std::cout << "e - 1 = " << acc_result << std::endl;
    double t_result = T(f, 0.0, 1.0, 1000);
    std::cout << "  T   = " << t_result << ", diff = " << t_result - acc_result << std::endl;
    double s_result = S(f, 0.0, 1.0, 8);
    std::cout << "  S   = " << s_result << ", diff = " << s_result - acc_result << std::endl;
    double l_result = L(f, 0.0, 1.0, eps);
    std::cout << "  L   = " << l_result << ", diff = " << l_result - acc_result << std::endl;
    return 0;
}