#include <iostream>
#define _USE_MATH_DEFINES
#include <cmath>

using namespace std;

double f(int m, double x, double t)
{
    return cos(m * t - x * sin(t));
}

double besselJ(int m, double x)
{
    int i;
    double temp = 0;
    double dt = M_PI / 100;
    for (i = 1; i <= 100; i++)
    {
        temp += (dt / 6) * (f(m, x, (i - 1) * dt) + 4 * f(m, x, (2 * i - 1) * dt / 2) + f(m, x, i * dt));
    }
    return temp / M_PI;
}

int main(int argc, const char *argv[])
{
    int i, n = 10000;
    double dx = 2 * M_PI / n;
    // double delta_x = 1e-5;
    double J0_der;

    for (i = 0; i < n; i++)
    {
        J0_der = (besselJ(0, i * dx + dx/2) - besselJ(0, i * dx - dx/2)) / (dx);
        cout << fabs(J0_der + besselJ(1, i * dx)) << "\n";
    }

    return 0;
}
