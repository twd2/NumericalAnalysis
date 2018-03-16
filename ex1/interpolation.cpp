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
        for (int i = 0; i < _samples->size(); ++i)
        {
            double divisor = 1;
            for (int j = 0; j < _samples->size(); ++j)
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
        for (int i = 0; i < _samples->size(); ++i)
        {
            double dividend = 1;
            for (int j = 0; j < _samples->size(); ++j)
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

class spline
{

};

int main()
{
    std::cout << std::fixed << std::setprecision(20);
    int n = 10;
    auto f = [](double x) { return (double)1 / (1 + 16 * x * x); };

    std::shared_ptr<samples_t> samples = std::make_shared<samples_t>(n + 1);
    for (int i = 0; i <= n; ++i)
    {
        double x = -5 + (double)i * 10 / n;
        (*samples)[i] = std::make_pair<double, double>((double)(x), f(x));
        std::cout << std::fixed << std::setprecision(2) << (*samples)[i].first << " "
                  << std::fixed << std::setprecision(20)
                  << (*samples)[i].second << std::endl;
    }
    
    lagrange la(samples);
    for (double x = -5.0; x <= 5.0; x += 0.0001)
    {
        std::cout << std::fixed << std::setprecision(2) << x << " "
                  << std::fixed << std::setprecision(20) << la(x) << " " << f(x) << std:: endl;
    }

    return 0;
}