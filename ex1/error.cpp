#include <cmath>
#include <iostream>
#include <iomanip>

int main()
{
    int n = 0, sign = -1;
    float xn = 0, eps = 0.5e-4, ln2 = 0.693147190546;
    while (std::abs(xn - ln2) > eps)
    {
        ++n;
        sign = -sign;
        xn += (decltype(xn))sign / (decltype(xn))n;
        std::cout << std::fixed << std::setprecision(20) << n << " " << xn << " " << (xn - ln2) << std::endl;
    }
    std::cout << "n: " << n << std::endl
              << "n expected: " << (int)(1 / eps - 1) << std::endl;
    // n > 1 / eps - 1
    return 0;
}