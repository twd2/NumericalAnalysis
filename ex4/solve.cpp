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
                long double sum = 0.0;
                for (int k = 0; k < r; ++k)
                {
                    sum += l(r, k) * u(k, j);
                }
                u(r, j) = a(r, j) - sum;
            }
            for (int i = r + 1; i < n; ++i)
            {
                long double sum = 0.0;
                for (int k = 0; k < r; ++k)
                {
                    sum += l(i, k) * u(k, r);
                }
                l(i, r) = (a(i, r) - sum) / u(r, r);
            }
        }
    }

    void sqrt(matrix &l, matrix &d) const
    {
        // 改进的平方根算法
        if (_row != _col)
        {
            exit(1);
        }

        const int n = _row;
        const matrix &a = *this;
        l = ident(n);
        d = zeros(n);

        matrix t = zeros(n);
        d(0, 0) = a(0, 0);
        for (int i = 1; i < n; ++i)
        {
            for (int j = 0; j < i; ++j)
            {
                long double sum = 0.0;
                for (int k = 0; k < j; ++k)
                {
                    sum += t(i, k) * l(j, k);
                }
                t(i, j) = a(i, j) - sum;
                l(i, j) = t(i, j) / d(j, j);
            }
            long double sum = 0.0;
            for (int k = 0; k < i; ++k)
            {
                sum += t(i, k) * l(i, k);
            }
            d(i, i) = a(i, i) - sum;
        }
    }

    matrix solve_as_l(const matrix &b) const
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
            long double sum = 0.0;
            for (int k = 0; k < i; ++k)
            {
                sum += l(i, k) * y._v[k];
            }
            y._v[i] = b._v[i] - sum;
        }
        return y;
    }

    matrix solve_as_u(const matrix &y) const
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
            long double sum = 0.0;
            for (int k = i + 1; k < n; ++k)
            {
                sum += u(i, k) * x._v[k];
            }
            x._v[i] = (y._v[i] - sum) / u(i, i);
        }
        return x;
    }

    matrix solve_as_lt_with_d(const matrix &y, const matrix &d) const
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
        const matrix &l = *this;
        matrix x(n, 1);
        x._v[n - 1] = y._v[n - 1] / d(n - 1, n - 1);
        for (int i = n - 2; i >= 0; --i)
        {
            long double sum = 0.0;
            for (int k = i + 1; k < n; ++k)
            {
                sum += l(k, i) * x._v[k];
            }
            x._v[i] = y._v[i] / d(i, i) - sum;
        }
        return x;
    }

    long double norm_inf() const
    {
        long double norm = _v[0];
        for (int i = 0; i < _row; ++i)
        {
            long double sum = 0.0;
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

    long double norm_row() const
    {
        return norm_inf();
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

    long double cond_inf() const
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

void cond_test()
{
    // (1)
    std::cout << std::scientific << std::setprecision(10);
    std::cout << "cond(H3)_inf = " << (double)matrix::hilbert(3).cond_inf() << std::endl;
    std::cout << "cond(H4)_inf = " << (double)matrix::hilbert(4).cond_inf() << std::endl;
}

void lu_test(const int n)
{
    // (2), (4), (5)
    std::cout << "========== LU ==========" << std::endl;
    matrix h = matrix::hilbert(n);
    std::cout << std::scientific << std::setprecision(10);
    std::cout << "cond(H" << n << ")_inf = " << (double)h.cond_inf() << std::endl;
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "H" << n << " = " << std::endl; h.print(); std::cout << std::endl;
    // LU分解
    matrix l, u;
    h.lu(l, u);
    std::cout << "L = " << std::endl; l.print(); std::cout << std::endl;
    std::cout << "U = " << std::endl; u.print(); std::cout << std::endl;
    // 验证
    std::cout << "L * U = " << std::endl; (l * u).print(); std::cout << std::endl;
    matrix x(n, 1); // 构造精确解
    for (int i = 0; i < n; ++i)
    {
        x(i, 0) = 1.0;
    }
    std::cout << std::fixed << std::setprecision(20);
    std::cout << "x = " << std::endl; x.print(); std::cout << std::endl;
    matrix b = h * x;
    std::cout << "b = H" << n << " * x = " << std::endl; b.print(); std::cout << std::endl;
    // 用LU分解结果计算近似解\hat x
    matrix x_hat = u.solve_as_u(l.solve_as_l(b));
    std::cout << "x_hat = " << std::endl; x_hat.print(); std::cout << std::endl;
    matrix r = b - h * x_hat, delta = x_hat - x;
    std::cout << std::scientific << std::setprecision(10);
    std::cout << "norm_inf(r) = " << (double)r.norm_inf() << std::endl;
    std::cout << "norm_inf(delta) = " << (double)delta.norm_inf() << std::endl;
    std::cout << "Applying 1e-7 jitter..." << std::endl;
    matrix b1 = b;
    b1(6, 0) += 1e-7;
    std::cout << std::fixed << std::setprecision(20);
    std::cout << "b1 = " << std::endl; b1.print(); std::cout << std::endl;
    matrix x1_hat = u.solve_as_u(l.solve_as_l(b1));
    std::cout << "x1_hat = " << std::endl; x1_hat.print(); std::cout << std::endl;
    matrix r1 = b1 - h * x1_hat, delta1 = x1_hat - x;
    std::cout << std::scientific << std::setprecision(10);
    std::cout << "norm_inf(r1) = " << (double)r1.norm_inf() << std::endl;
    std::cout << "norm_inf(delta1) = " << (double)delta1.norm_inf() << std::endl;
}

void sqrt_test(const int n)
{
    // (3)
    std::cout << "========== Cholesky ==========" << std::endl;
    matrix h = matrix::hilbert(n);
    std::cout << std::fixed << std::setprecision(6);
    // Cholesky分解
    matrix l, d;
    h.sqrt(l, d);
    std::cout << "L = " << std::endl; l.print(); std::cout << std::endl;
    std::cout << "D = " << std::endl; d.print(); std::cout << std::endl;
    // 验证
    std::cout << "L * D * LT = " << std::endl; (l * d * l.transpose()).print(); std::cout << std::endl;
    matrix x(n, 1); // 构造精确解
    for (int i = 0; i < n; ++i)
    {
        x(i, 0) = 1.0;
    }
    std::cout << std::fixed << std::setprecision(20);
    std::cout << "x = " << std::endl; x.print(); std::cout << std::endl;
    matrix b = h * x;
    std::cout << "b = H" << n << " * x = " << std::endl; b.print(); std::cout << std::endl;
    // 用Cholesky分解结果计算近似解\hat x
    matrix x_hat = l.solve_as_lt_with_d(l.solve_as_l(b), d);
    std::cout << "x_hat = " << std::endl; x_hat.print(); std::cout << std::endl;
    matrix r = b - h * x_hat, delta = x_hat - x;
    std::cout << std::scientific << std::setprecision(10);
    std::cout << "norm_inf(r) = " << (double)r.norm_inf() << std::endl;
    std::cout << "norm_inf(delta) = " << (double)delta.norm_inf() << std::endl;
}

void benchmark(const int n, const int iter = 10000)
{
    using namespace std::chrono;
    std::cout << "========== benchmark ==========" << std::endl;
    matrix h = matrix::hilbert(n);
    matrix x(n, 1); // 构造精确解
    for (int i = 0; i < n; ++i)
    {
        x(i, 0) = 1.0;
    }
    matrix b = h * x;

	steady_clock::time_point start_time = steady_clock::now();
    steady_clock::time_point end_time = steady_clock::now();
    double bias = duration_cast<duration<double> >(end_time - start_time).count();

    double lu_time;
    start_time = steady_clock::now();
    for (int i = 0; i < iter; ++i)
    {
        matrix l, u;
        h.lu(l, u);
        matrix x_hat = u.solve_as_u(l.solve_as_l(b));
    }
    end_time = steady_clock::now();
    lu_time = duration_cast<duration<double> >(end_time - start_time).count() - bias;
    lu_time /= iter;

    double sqrt_time;
    start_time = steady_clock::now();
    for (int i = 0; i < iter; ++i)
    {
        matrix l, d;
        h.sqrt(l, d);
        matrix x_hat = l.solve_as_lt_with_d(l.solve_as_l(b), d);
    }
    end_time = steady_clock::now();
    sqrt_time = duration_cast<duration<double> >(end_time - start_time).count() - bias;
    sqrt_time /= iter;

    std::cout << "n = " << n << std::endl;
    std::cout << "iter = " << iter << std::endl;
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "LU time = " << lu_time * 1e9 << " ns" << std::endl;
    std::cout << "Cholesky time = " << sqrt_time * 1e9 << " ns" << std::endl;
}

void n_test(const int n)
{
    std::cout << std::scientific << std::setprecision(10);
    std::cout << "n = " << n << std::endl;
    matrix h = matrix::hilbert(n);
    // LU分解
    matrix l, u;
    h.lu(l, u);
    matrix x(n, 1); // 构造精确解
    for (int i = 0; i < n; ++i)
    {
        x(i, 0) = 1.0;
    }
    matrix b = h * x;
    // 用LU分解结果计算近似解\hat x
    matrix x_hat = u.solve_as_u(l.solve_as_l(b));
    matrix r = b - h * x_hat, delta = x_hat - x;
    std::cout << "norm_inf(r) = " << (double)r.norm_inf() << std::endl;
    std::cout << "norm_inf(delta) = " << (double)delta.norm_inf() << std::endl;
}

int main()
{
    cond_test();
    lu_test(10);
    sqrt_test(10);
    benchmark(100, 10000);
    std::cout << "========== n ==========" << std::endl;
    for (int n = 1; n <= 100; ++n)
    {
        n_test(n);
    }
    return 0;
}