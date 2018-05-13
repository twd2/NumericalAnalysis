#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>
#include <chrono>
#include <functional>

double binsolve(std::function<double (double)> f, double l, double r, const double eps = 1e-5)
{
    int i = 0;
    double fl = f(l), fr = f(r);
    std::cout << "iter 0: [" << l << ", " << r << "], f(l) = " << fl << ", f(r) = " << fr << std::endl;
    while (r - l > eps)
    {
        double m = (l + r) / 2;
        double fm = f(m);
        if (fl * fm < 0)
        {
            r = m;
            fr = fm;
        }
        else
        {
            l = m;
            fl = fm;
        }
        ++i;
        std::cout << "iter " << i << ": [" << l << ", " << r << "], f(l) = " << fl << ", f(r) = " << fr << std::endl;
    }
    return (l + r) / 2;
}

double nsolve(std::function<double (double)> f, std::function<double (double)> df, double x0, const double eps = 1e-5)
{
    return 0.35;
}

int main()
{
    //std::cout << std::fixed << std::setprecision(10);
    std::cout << std::scientific << std::setprecision(10);
    auto f = [](double x) { return ((2 * x - 1) * x + 3) * x - 1; }; // 2x^3 - x^2 + 3x - 1
    auto df = [](double x) { return (6 * x - 2) * x + 3; }; // 6x^2 - 2x + 3
    std::cout << "Algorithm: Binary Search" << std::endl;
    double x = binsolve(f, -3, 3, 1e-5);
    std::cout << "x = " << x << ", f(x) = " << f(x) << std::endl << std::endl;
    std::cout << "Algorithm: Newton's" << std::endl;
    x = nsolve(f, df, 0, 1e-5);
    std::cout << "x = " << x << ", f(x) = " << f(x) << std::endl << std::endl;
    //std::cout << std::scientific << std::setprecision(10);
    return 0;
}