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

    }

    poly operator*(const poly &b) const
    {
        poly result(this->n() * b.n() + 1);
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

    double length2(const std::vector<double> &x) const
    {
        double result = 0;
        for (std::size_t i = 0; i < x.size(); ++i)
        {
            double y = (*this)(x[i]);
            result += y * y;
        }
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
    std::vector<double> x { -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0 },
                        y { -4.467, -0.452, 0.551, 0.048, -0.447, 0.549, 4.552 };
    /*poly p1({-1, 1}), p2({1, 1});
    p1.print(); std::cout << std::endl;
    //p2 = p1;
    p2 = p1 * p2;
    p2.print(); std::cout << std::endl;
    (poly({-1,0}) * (p1 + p2)).print(); std::cout << std::endl;
    std::cout << p2(3.3) << " " << p2(2) << std::endl;*/
    return 0;
}