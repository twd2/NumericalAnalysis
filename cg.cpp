#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>
#include <chrono>

class matrix
{
private:
    std::vector<long double> _v;
    int _row, _col;
public:
    matrix(const matrix &m)
        : _v(m._v), _row(m._row), _col(m._col)
    {

    }

    matrix(matrix &&m)
        : _v(std::move(m._v)), _row(m._row), _col(m._col)
    {

    }

    matrix &operator=(const matrix &m)
    {
        _v = m._v;
        _row = m._row;
        _col = m._col;
        return *this;
    }

    matrix &operator=(matrix &&m)
    {
        _v = std::move(m._v);
        _row = m._row;
        _col = m._col;
        return *this;
    }

    matrix(int row, int col)
        : _v(row * col, 0), _row(row), _col(col)
    {

    }

    matrix()
        : matrix(0, 0)
    {

    }

    matrix(int row, int col, const std::vector<long double> &v)
        : _v(v), _row(row), _col(col)
    {

    }

    matrix(int row, int col, std::vector<long double> &&v)
        : _v(std::move(v)), _row(row), _col(col)
    {

    }

    void print()
    {
        int k = 0;
        for (int i = 0; i < _row; ++i)
        {
            for (int j = 0; j < _col; ++j)
            {
                std::cout << (double)_v[k] << " ";
                ++k;
            }
            std::cout << std::endl;
        }
    }

    long double &operator()(int i, int j)
    {
        return _v[i * _col + j];
    }

    long double operator()(int i, int j) const
    {
        return _v[i * _col + j];
    }

    int row() const
    {
        return _row;
    }

    int col() const
    {
        return _col;
    }

    matrix operator*(const matrix &b) const
    {
        if (_col != b._row)
        {
            exit(1);
        }

        const matrix &a = *this;

        matrix c(_row, b._col);
        for (int i = 0; i < c._row; ++i)
        {
            for (int j = 0; j < c._col; ++j)
            {
                long double sum = 0.0;
                for (int k = 0; k < _col; ++k)
                {
                    sum += a(i, k) * b(k, j);
                }
                c(i, j) = sum;
            }
        }
        return c;
    }

    matrix operator+(const matrix &b) const
    {
        if (_row != b._row || _col != b._col)
        {
            exit(1);
        }

        const matrix &a = *this;

        matrix c(_row, _col);
        for (int i = 0; i < _row; ++i)
        {
            for (int j = 0; j < _col; ++j)
            {
                c(i, j) = a(i, j) + b(i, j);
            }
        }
        return c;
    }

    matrix operator-(const matrix &b) const
    {
        if (_row != b._row || _col != b._col)
        {
            exit(1);
        }

        const matrix &a = *this;

        matrix c(_row, _col);
        for (int i = 0; i < _row; ++i)
        {
            for (int j = 0; j < _col; ++j)
            {
                c(i, j) = a(i, j) - b(i, j);
            }
        }
        return c;
    }

    matrix operator*(const double b) const
    {
        const matrix &a = *this;

        matrix c(_row, _col);
        for (int i = 0; i < _row; ++i)
        {
            for (int j = 0; j < _col; ++j)
            {
                c(i, j) = b * a(i, j);
            }
        }
        return c;
    }

    matrix transpose() const
    {
        matrix result(_col, _row);
        for (int i = 0; i < _row; ++i)
        {
            for (int j = 0; j < _col; ++j)
            {
                result(j, i) = (*this)(i, j);
            }
        }
        return result;
    }

    double as_double() const
    {
        if (_row != 1 || _col != 1)
        {
            exit(1);
        }
        return _v[0];
    }

    static matrix zeros(int n);
    static matrix ident(int n);
};

matrix matrix::zeros(int n)
{
    return matrix(n, n);
}

matrix matrix::ident(int n)
{
    matrix m(n, n);
    for (int i = 0; i < n; ++i)
    {
        m(i, i) = 1.0;
    }
    return m;
}

int main()
{
    std::cout << std::scientific << std::setprecision(10);
    matrix A(3, 3, { 4, 3, 0, 3, 4, -1, 0, -1, 4 }),
           y(3, 1, { 3, 5, -5 });
    matrix x(3, 1);
    matrix r0, r, p;
    double a, b;
    std::cout << "x = " << std::endl; x.print(); std::cout << std::endl;
    r = y - A * x;
    std::cout << "r = " << std::endl; r.print(); std::cout << std::endl;
    p = r;
    std::cout << "p = " << std::endl; p.print(); std::cout << std::endl;
    for (int i = 0; i < 5; ++i)
    {
        std::cout << (r.transpose() * r).as_double() << std::endl;
        std::cout << (p.transpose() * A * p).as_double() << std::endl;
        a = (r.transpose() * r).as_double() / (p.transpose() * A * p).as_double();
        std::cout << "a = " << a << std::endl << std::endl;
        x = x + p * a;
        std::cout << "x = " << std::endl; x.print(); std::cout << std::endl;
        r0 = r;
        r = y - A * x;
        std::cout << "r = " << std::endl; r.print(); std::cout << std::endl;
        b = (r.transpose() * r).as_double() / (r0.transpose() * r0).as_double();
        std::cout << "b = " << b << std::endl << std::endl;
        p = r + p * b;
        std::cout << "p = " << std::endl; p.print(); std::cout << std::endl;
    }
    return 0;
}