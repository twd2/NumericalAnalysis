#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>

#define DEBUG 1

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

    void print() const
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

    matrix eigenvec(const double eps = 1e-5) const
    {
        if (_row != _col)
        {
            exit(1);
        }

        const int n = _row;
        matrix vec(n, 1);
        for (int i = 0; i < n; ++i)
        {
            vec._v[i] = 1.0; // FIXME: this is an arbitrary value.
        }
#ifdef DEBUG
        int iter = 0;
#endif
        double lambda, lambda1 = vec._v[vec.absmax_index()];
        do
        {
            vec = vec * (1 / lambda1);
#ifdef DEBUG
            std::cout << "iter = " << iter << ", " << "vec = " << std::endl; vec.print(); std::cout << std::endl;
#endif
            vec = (*this) * vec;
            lambda = lambda1;
            lambda1 = vec._v[vec.absmax_index()];
#ifdef DEBUG
            ++iter;
#endif
        } while (std::abs(lambda1 - lambda) >= eps);
        return vec * (1 / lambda1);
    }

    int absmax_index() const
    {
        double maxval = std::abs(_v[0]);
        int maxid = 0;
        for (int i = 1; i < _row; ++i)
        {
            double av = std::abs(_v[i]);
            if (av > maxval)
            {
                maxval = av;
                maxid = i;
            }
        }
        return maxid;
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

void do_eigen(const matrix &m)
{
    std::cout << "m = " << std::endl; m.print(); std::cout << std::endl;
    matrix ev = m.eigenvec();
    std::cout << "ev = " << std::endl; ev.print(); std::cout << std::endl;
    matrix mev = m * ev;
    double lambda = mev(mev.absmax_index(), 0);
    std::cout << "lambda = " << lambda << std::endl << std::endl;
}

int main()
{
    std::cout << std::fixed << std::setprecision(10);
    do_eigen(matrix(3, 3, { 5, -4, 1, -4, 6, -4, 1, -4, 7 }));
    do_eigen(matrix(4, 4, { 25, -41, 10, -6, -41, 68, -17, 10, 10, -17, 5, -3, -6, 10, -3, 2 }));
    return 0;
}