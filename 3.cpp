#include <iostream>
#define _USE_MATH_DEFINES
#include <cmath>

const int INTERV = 1<<1;

using namespace std;

double func(const double x, const int which)
{
    if (which == 1)
    {
        return 1 / (1 + x * x); // first function
    }
    else if (which == 2)
    {
        return pow(x, (1. / 3)) * exp(sin(x)); // second function
    }
}

double trapezoidal(const int which, const int intervals, const double left, const double right)
{
    double temp = 0;
    double dx = (right - left) / intervals;
    int i;

    for (i = 0; i < intervals; i++)
    {
        temp += .5 * dx * (func(left + dx * i, which) + func(left + dx * (i + 1), which));
    }
    return temp;
}

double simpson(const int which, const int intervals, const double left, const double right)
{
    double temp = 0;
    double dx = (right - left) / intervals;
    int i;
    double a, b;

    for (i = 1; i <= intervals; i++)
    {
        a = left + (i - 1) * dx;
        b = left + i * dx;
        temp += (dx / 6) * (func(a, which) + 4 * func((a + b) / 2, which) + func(b, which));
    }

    return temp;
}

int main() {
    cout << "Exact value | Trapezoidal | Simpson\n" << endl;
    cout << "first func:\n" << M_PI_2 << "   " << trapezoidal(1, INTERV, -1, 1) << "   " << simpson(1, INTERV, -1, 1) << endl;
    cout << "second func:\n" << "1.29587" << "   " << trapezoidal(2, INTERV, 0, 1) << "   " << simpson(2, INTERV, 0, 1) << endl;
    return 0;
}
