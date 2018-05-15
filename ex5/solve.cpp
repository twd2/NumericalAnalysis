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
    std::cout << "iter 0: x in [" << l << ", " << r << "], f(l) = " << fl << ", f(r) = " << fr << std::endl;
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
        std::cout << "iter " << i << ": x in [" << l << ", " << r << "], f(l) = " << fl << ", f(r) = " << fr << std::endl;
    }
    return (l + r) / 2;
}

double nsolve(std::function<double (double)> f, std::function<double (double)> df, double x0, const double eps = 1e-5)
{
    int i = 0;
    double x = x0;
    double fx = f(x), dfx = df(x);
    double delta = fx / dfx;
    std::cout << "iter 0: x = " << x << ", f(x) = " << fx << ", f'(x) = " << dfx << std::endl;
    while (std::abs(delta) >= eps)
    {
        x -= delta;
        fx = f(x);
        dfx = df(x);
        delta = fx / dfx;
        ++i;
        std::cout << "iter " << i << ": x = " << x << ", f(x) = " << fx << ", f'(x) = " << dfx << std::endl;
    }
    return x;
}

double itersolve(std::function<double (double)> phi, double x0, const int iter = 10)
{
    double x = x0;
    for (int i = 0; i <= iter; ++i)
    {
        std::cout << "iter " << i << ": x = " << x << std::endl;
        x = phi(x);
    }
    return x;
}

int main()
{
    std::cout << std::scientific << std::setprecision(10);
    {
        std::cout << "========== Problem 1 ==========" << std::endl;
        auto f = [](double x) { return ((2 * x - 1) * x + 3) * x - 1; }; // 2x^3 - x^2 + 3x - 1
        auto df = [](double x) { return (6 * x - 2) * x + 3; }; // 6x^2 - 2x + 3
        std::cout << "Algorithm: Binary Search" << std::endl;
        double x = binsolve(f, -3, 3, 1e-5);
        std::cout << "x = " << x << ", f(x) = " << f(x) << std::endl << std::endl;
        std::cout << "Algorithm: Newton's" << std::endl;
        x = nsolve(f, df, 0, 1e-5);
        std::cout << "x = " << x << ", f(x) = " << f(x) << std::endl << std::endl;
    }

    {
        std::cout << "========== Problem 2.1 ==========" << std::endl;
        auto f = [](double x) { return (2 * x * x - 1) * x - 1; }; // 2x^3 - x - 1
        auto phi1 = [](double x) { return std::cbrt((x + 1) / 2); };
        double x = itersolve(phi1, 0);
        std::cout << "x = " << x << ", f(x) = " << f(x) << std::endl << std::endl;
        auto phi2 = [](double x) { return 2 * x * x * x - 1; };
        x = itersolve(phi2, 0);
        std::cout << "x = " << x << ", f(x) = " << f(x) << std::endl << std::endl;
    }

    {
        std::cout << "========== Problem 2.2 ==========" << std::endl;
        auto f = [](double x) { return (x * x - 1) * x - 1; }; // x^3 - x - 1
        auto df = [](double x) { return 3 * x * x - 1; }; // 3x^2 - 1
        double x;
        std::cout << "Newton's x0 = 1.5" << std::endl;
        x = nsolve(f, df, 1.5, 1e-5);
        std::cout << "x = " << x << ", f(x) = " << f(x) << std::endl << std::endl;
        std::cout << "Newton's x0 = 0.0" << std::endl;
        x = nsolve(f, df, 0, 1e-5);
        std::cout << "x = " << x << ", f(x) = " << f(x) << std::endl << std::endl;
    }
    return 0;
}