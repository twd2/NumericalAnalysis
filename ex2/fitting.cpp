#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>

class poly
{
private:
    std::vector<double> _coeff;
public:
    poly()
        : _coeff()
    {

    }

    explicit poly(int n)
        : _coeff(n + 1, 0.0)
    {

    }

    explicit poly(const std::vector<double> &coeff)
        : _coeff(coeff)
    {

    }

    poly(const poly &b)
        : _coeff(b._coeff)
    {

    }

    poly &operator=(const poly &b)
    {
        _coeff = b._coeff;
        return *this;
    }

    void normalize()
    {
        while (_coeff.size() >= 2 && _coeff.back() == 0.0)
        {
            _coeff.pop_back();
        }
    }

    poly operator*(const poly &b) const
    {
        poly result(this->n() + b.n() + 1);
        for (int i = 0; i <= this->n(); ++i)
        {
            for (int j = 0; j <= b.n(); ++j)
            {
                result._coeff[i + j] += _coeff[i] * b._coeff[j];
            }
        }
        result.normalize();
        return result;
    }

    poly operator+(const poly &b) const
    {
        poly result(std::max(this->n(), b.n()) + 1);
        for (int i = 0; i <= result.n(); ++i)
        {
            result._coeff[i] = (*this)[i] + b[i];
        }
        result.normalize();
        return result;
    }

    double dot(const std::vector<double> &y, const std::vector<double> &x) const
    {
        double result = 0;
        for (std::size_t i = 0; i < x.size(); ++i)
        {
            result += (*this)(x[i]) * y[i];
        }
        return result;
    }

    double dot(const poly &b, const std::vector<double> &x) const
    {
        double result = 0;
        for (std::size_t i = 0; i < x.size(); ++i)
        {
            result += (*this)(x[i]) * b(x[i]);
        }
        return result;
    }

private:
    double _length2_cache = 0.0;
    const std::vector<double> *_last_x = nullptr;
public:
    double length2(const std::vector<double> &x)
    {
        if (&x == _last_x)
        {
            return _length2_cache;
        }
        double result = 0;
        for (std::size_t i = 0; i < x.size(); ++i)
        {
            double y = (*this)(x[i]);
            result += y * y;
        }
        _length2_cache = result;
        _last_x = &x;
        return result;
    }

    double xdot(const std::vector<double> &x) const
    {
        // xP(x) dot P(x)
        double result = 0;
        for (std::size_t i = 0; i < x.size(); ++i)
        {
            double y = (*this)(x[i]);
            result += x[i] * y * y;
        }
        return result;
    }

    double operator()(double x) const
    {
        // 秦九韶算法
        double result = _coeff[n()];
        for (int i = n() - 1; i >= 0; --i)
        {
            result = result * x + _coeff[i];
        }
        return result;
    }

    int n() const
    {
        return _coeff.size() - 1;
    }

    double operator[](std::size_t i) const
    {
        if (i >= _coeff.size())
        {
            return 0;
        }
        return _coeff[i];
    }

    double &operator[](std::size_t i)
    {
        return _coeff[i];
    }

    void print() const
    {
        bool printed = false;
        for (int i = n(); i >= 0; --i)
        {
            bool sign = _coeff[i] < 0;
            if (_coeff[i] == 0)
            {
                continue;
            }
            if (printed && !sign)
            {
                std::cout << " + ";
            }
            if (sign)
            {
                std::cout << " - ";
            }
            if ((_coeff[i] != 1.0 && _coeff[i] != -1.0) || i == 0)
            {
                std::cout << std::abs(_coeff[i]);
            }
            
            if (i == 1)
            {
                std::cout << "x";
            }
            else if (i == 0)
            {

            }
            else
            {
                std::cout << "x^" << i;
            }
            printed = true;
        }
        if (!printed)
        {
            std::cout << "0";
        }
    }
};

int main()
{
    std::cout << std::fixed << std::setprecision(20);
    const int n = 3;
    std::vector<double> x { -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0 },
                        y { -4.467, -0.452, 0.551, 0.048, -0.447, 0.549, 4.552 };

    // 计算正交基
    double alpha, beta;
    poly p[n + 1];
    p[0] = poly(std::vector<double> { 1.0 });
    alpha = p[0].xdot(x) / x.size();
    p[1] = poly({ -alpha, 1.0 });
    for (int i = 2; i <= n; ++i)
    {
        alpha = p[i - 1].xdot(x) / p[i - 1].length2(x);
        beta = p[i - 1].length2(x) / p[i - 2].length2(x);
        p[i] = poly({ -alpha, 1.0 }) * p[i - 1] + poly(std::vector<double> { -beta }) * p[i - 2];
    }

    std::cout << "Bases:" << std::endl;
    for (int i = 0; i <= n; ++i)
    {
        p[i].print();
        std::cout << std::endl;
    }
    std::cout << std::endl;

    poly fit;
    for (int i = 0; i <= n; ++i)
    {
        std::cout << "n = " << i << std::endl;
        double coeff = p[i].dot(y, x) / p[i].length2(x);
        fit = fit + poly(std::vector<double> { coeff }) * p[i];
        std::cout << "pn = "; fit.print(); std::cout << std::endl;
        std::cout << "x\ty\tpn(x)\tdiff" << std::endl;
        for (std::size_t j = 0; j < x.size(); ++j)
        {
            std::cout << x[j] << "\t" << y[j] << "\t"
                      << fit(x[j]) << "\t" << fit(x[j]) - y[j] << std::endl;
        }
        std::cout << std::endl;
    }
    return 0;
}