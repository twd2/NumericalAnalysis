#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>

class matrix
{
private:
    std::vector<double> _v;
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

    matrix(int row, int col, const std::vector<double> &v)
        : _v(v), _row(row), _col(col)
    {

    }

    matrix(int row, int col, std::vector<double> &&v)
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
                std::cout << _v[k] << " ";
                ++k;
            }
            std::cout << std::endl;
        }
    }

    double &operator()(int i, int j)
    {
        return _v[i * _col + j];
    }

    double operator()(int i, int j) const
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
                double sum = 0.0;
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

    void lu(matrix &l, matrix &u) const
    {
        if (_row != _col)
        {
            exit(1);
        }

        const int n = _row;
        const matrix &a = *this;

        l = ident(n);
        u = zeros(n);
        for (int j = 0; j < n; ++j)
        {
            u(0, j) = a(0, j);
        }
        for (int i = 1; i < n; ++i)
        {
            l(i, 0) = a(i, 0) / u(0, 0);
        }
        for (int r = 1; r < n; ++r)
        {
            for (int j = r; j < n; ++j)
            {
                double sum = 0.0;
                for (int k = 0; k < r; ++k)
                {
                    sum += l(r, k) * u(k, j);
                }
                u(r, j) = a(r, j) - sum;
            }
            for (int i = r + 1; i < n; ++i)
            {
                double sum = 0.0;
                for (int k = 0; k < r; ++k)
                {
                    sum += l(i, k) * u(k, r);
                }
                l(i, r) = (a(i, r) - sum) / u(r, r);
            }
        }
    }

    matrix solve_as_l(const matrix &b)
    {
        if (_row != _col)
        {
            exit(1);
        }

        if (b._row != _row)
        {
            exit(1);
        }

        if (b._col != 1)
        {
            exit(1);
        }

        const int n = _row;
        const matrix &l = *this;
        matrix y(n, 1);
        y._v[0] = b._v[0];
        for (int i = 1; i < n; ++i)
        {
            double sum = 0.0;
            for (int k = 0; k < i; ++k)
            {
                sum += l(i, k) * y._v[k];
            }
            y._v[i] = b._v[i] - sum;
        }
        return y;
    }

    matrix solve_as_u(const matrix &y)
    {
        if (_row != _col)
        {
            exit(1);
        }

        if (y._row != _row)
        {
            exit(1);
        }

        if (y._col != 1)
        {
            exit(1);
        }

        const int n = _row;
        const matrix &u = *this;
        matrix x(n, 1);
        x._v[n - 1] = y._v[n - 1] / u(n - 1, n - 1);
        for (int i = n - 2; i >= 0; --i)
        {
            double sum = 0.0;
            for (int k = i + 1; k < n; ++k)
            {
                sum += u(i, k) * x._v[k];
            }
            x._v[i] = (y._v[i] - sum) / u(i, i);
        }
        return x;
    }

    double norm_inf() const
    {
        double norm = _v[0];
        for (int i = 0; i < _row; ++i)
        {
            double sum = 0.0;
            for (int j = 0; j < _col; ++j)
            {
                sum += std::abs((*this)(i, j));
            }
            if (sum > norm)
            {
                norm = sum;
            }
        }
        return norm;
    }

    double norm_row() const
    {
        return norm_inf();
    }

    matrix inverse() const
    {
        if (_row != _col)
        {
            exit(1);
        }

        const int n = _row;
        matrix l, u;
        lu(l, u);
        matrix inv(n, n);
        matrix b(n, 1);
        for (int j = 0; j < n; ++j)
        {
            for (int i = 0; i < n; ++i)
            {
                if (i == j)
                {
                    b._v[i] = 1.0;
                }
                else
                {
                    b._v[i] = 0.0;
                }
            }
            matrix x = u.solve_as_u(l.solve_as_l(b));
            for (int i = 0; i < n; ++i)
            {
                inv(i, j) = x._v[i];
            }
        }
        return inv;
    }

    double cond_inf() const
    {
        matrix inv = inverse();
        return inv.norm_inf() * norm_inf();
    }

    static matrix hilbert(int n);
    static matrix zeros(int n);
    static matrix ident(int n);
};

matrix matrix::hilbert(int n)
{
    matrix m(n, n);
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            m(i, j) = 1.0 / (i + j + 1);
        }
    }
    return m;
}

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
    // (1)
    std::cout << std::fixed << std::setprecision(20);
    std::cout << "cond(H3) = " << matrix::hilbert(3).cond_inf() << std::endl;
    std::cout << "cond(H4) = " << matrix::hilbert(4).cond_inf() << std::endl;
    // (2)
    const int n = 100;
    std::cout << "cond(H" << n << ") = " << matrix::hilbert(n).cond_inf() << std::endl;
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "H" << n << " = " << std::endl;
    matrix h = matrix::hilbert(n);
    h.print();
    std::cout << std::endl;
    matrix l, u;
    h.lu(l, u);
    std::cout << "L = " << std::endl;
    l.print();
    std::cout << std::endl;
    std::cout << "U = " << std::endl;
    u.print();
    std::cout << std::endl;
    std::cout << "L * U = " << std::endl;
    (l * u).print();
    std::cout << std::endl;
    matrix x(n, 1); // 精确解
    for (int i = 0; i < n; ++i)
    {
        x(i, 0) = 1.0;
    }
    std::cout << "x = " << std::endl;
    x.print();
    std::cout << std::endl;
    matrix b = h * x;
    std::cout << "b = h * x = " << std::endl;
    b.print();
    std::cout << std::endl;
    matrix x_hat = u.solve_as_u(l.solve_as_l(b));
    std::cout << "x_hat = " << std::endl;
    x_hat.print();
    std::cout << std::endl;
    matrix r = b - h * x_hat,
           delta = x_hat - x;
    std::cout << std::fixed << std::setprecision(20);
    std::cout << "norm_inf(r) = " << r.norm_inf() << std::endl;
    std::cout << "norm_inf(delta) = " << delta.norm_inf() << std::endl;
    return 0;
}