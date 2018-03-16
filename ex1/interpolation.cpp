#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>
#include <memory>

typedef std::vector<std::pair<double, double> > samples_t;

class lagrange
{
private:
    std::shared_ptr<samples_t> _samples;
    std::vector<double> _divisors;
public:
    lagrange(const std::shared_ptr<samples_t> &samples)
        : _samples(samples)
    {
        _divisors.resize(_samples->size());
        for (std::size_t i = 0; i < _samples->size(); ++i)
        {
            double divisor = 1;
            for (std::size_t j = 0; j < _samples->size(); ++j)
            {
                if (j == i)
                {
                    continue;
                }
                divisor *= (*_samples)[i].first - (*_samples)[j].first; // xk - xj
            }
            _divisors[i] = divisor;
        }
    }

    double operator()(double x) const
    {
        double result = 0;
        for (std::size_t i = 0; i < _samples->size(); ++i)
        {
            double dividend = 1;
            for (std::size_t j = 0; j < _samples->size(); ++j)
            {
                if (j == i)
                {
                    continue;
                }
                dividend *= x - (*_samples)[j].first; // x - xj
            }
            result += (*_samples)[i].second * dividend / _divisors[i];
        }
        return result;
    }
};

/*  追赶法解形如这样的方程：
    +-   -+    +- -+
    |bc   |    | f |
    |abc  |    | f |
    | abc |x = | f |
    |  abc|    | f |
    |   ab|    | f |
    +-   -+    +- -+
*/
std::vector<double> solve_equ(const std::vector<double> &a,
                              const std::vector<double> &b,
                              const std::vector<double> &c,
                              const std::vector<double> &f)
{
    const int n = b.size();
    std::vector<double> beta(n);

    beta[0] = c[0] / b[0];
    for (int i = 1; i < n - 1; ++i)
    {
        beta[i] = c[i] / (b[i] - a[i] * beta[i - 1]);
    }
    
    std::vector<double> xy(n);
    std::vector<double> &x = xy, &y = xy;
    
    y[0] = f[0] / b[0];
    for (int i = 1; i < n; ++i)
    {
        y[i] = (f[i] - a[i] * y[i - 1]) / (b[i] - a[i] * beta[i - 1]);
    }
    
    x[n - 1] = y[n - 1];
    for (int i = n - 2; i >= 0; --i)
    {
        x[i] = y[i] - beta[i] * x[i + 1];
    }
    return xy;
}

class spline
{
private:
    std::shared_ptr<samples_t> _samples;
    std::vector<double> _M;
    double _df0, _dfn;
public:
    spline(const std::shared_ptr<samples_t> &samples, double df0, double dfn)
        : _samples(samples), _M(samples->size()), _df0(df0), _dfn(dfn)
    {
        solve_m();
    }

    double operator()(double x) const
    {
        if (x == (*_samples)[_samples->size() - 1].first)
        {
            return (*_samples)[_samples->size() - 1].second;
        }
        for (std::size_t i = 0; i < _samples->size() - 1; ++i)
        {
            double x0 = (*_samples)[i].first, y0 = (*_samples)[i].second,
                   x1 = (*_samples)[i + 1].first, y1 = (*_samples)[i + 1].second;
            if (x0 <= x && x < x1)
            {
                double h = x1 - x0;
                // 书上式 (6.8)
                return _M[i] * (x1 - x)  * (x1 - x) * (x1 - x) / (6 * h) +
                       _M[i + 1] * (x - x0)  * (x - x0) * (x - x0) / (6 * h) +
                       (y0 - _M[i] * h * h / 6) * (x1 - x) / h +
                       (y1 - _M[i + 1] * h * h / 6) * (x - x0) / h;
            }
        }
        // ???
        return 0.0 / 0.0;
    }

private:
    void solve_m()
    {
        const int n = _samples->size() - 1;
        std::vector<double> a(_samples->size()), b(_samples->size(), 2.0),
                            c(_samples->size()), f(_samples->size());
        for (int i = 1; i <= n - 1; ++i)
        {
            double x0 = (*_samples)[i - 1].first, y0 = (*_samples)[i - 1].second,
                   x1 = (*_samples)[i + 0].first, y1 = (*_samples)[i + 0].second,
                   x2 = (*_samples)[i + 1].first, y2 = (*_samples)[i + 1].second;
            double h0 = x1 - x0, h1 = x2 - x1;
            // 书上式 (6.11)
            a[i] = h0 / (h0 + h1);
            c[i] = h1 / (h0 + h1);
            f[i] = 6.0 * ((y2 - y1) / (x2 - x1) - (y1 - y0) / (x1 - x0)) / (h0 + h1);
        }
        // 第一种边界
        c[0] = 1.0;
        // f[0] = 6.0 * ((y1 - y0) / (x1 - x0) - df0) / (x1 - x0);
        f[0] = 6.0 * (((*_samples)[1].second - (*_samples)[0].second) / ((*_samples)[1].first - (*_samples)[0].first) - _df0) / ((*_samples)[1].first - (*_samples)[0].first);
        a[n] = 1.0;
        // f[n] = 6.0 * (dfn - (yn - yn_1) / (xn - xn_1)) / (xn - xn_1);
        f[n] = 6.0 * (_dfn - ((*_samples)[n].second - (*_samples)[n - 1].second) / ((*_samples)[n].first - (*_samples)[n - 1].first)) / ((*_samples)[n].first - (*_samples)[n - 1].first);

        _M = solve_equ(a, b, c, f);
    }
};

int main()
{
    std::cout << std::fixed << std::setprecision(20);
    int n = 20;
    auto f = [](double x) { return (double)1 / (1 + 16 * x * x); };
    auto df = [](double x) { return (double)-(32 * x) / ((1 + 16 * x * x) * (1 + 16 * x * x)); };

    std::shared_ptr<samples_t> samples = std::make_shared<samples_t>(n + 1);
    for (int i = 0; i <= n; ++i)
    {
        double x = -5 + (double)i * 10 / n;
        (*samples)[i] = std::make_pair<double, double>((double)(x), f(x));
        std::cout << std::fixed << std::setprecision(2) << (*samples)[i].first << "\t"
                  << std::fixed << std::setprecision(20)
                  << (*samples)[i].second << std::endl;
    }
    
    lagrange la(samples);
    for (double x = -5.0; x <= 5.0; x += 0.01)
    {
        std::cout << std::fixed << std::setprecision(2) << x << "\t"
                  << std::fixed << std::setprecision(20) << la(x) << "\t" << f(x) << std:: endl;
    }

    std::cout << "--------" << std::endl;

    spline sp(samples, df((*samples)[0].first), df((*samples)[n].first));
    for (double x = -5.0; x <= 5.0; x += 0.01)
    {
        std::cout << std::fixed << std::setprecision(2) << x << "\t"
                  << std::fixed << std::setprecision(20) << sp(x) << "\t" << f(x) << std:: endl;
    }

    return 0;
}