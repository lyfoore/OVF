#include <iostream>
#define _USE_MATH_DEFINES
#include <cmath>

#define A 1
#define U0 1
#define epsilon 2e-15
#define ALPHA 1 // A*sqrt(U0)
#define E1 1 - M_PI *M_PI / (ALPHA * ALPHA)

using namespace std;

int n1 = 0;
int n2 = 0;
int n3 = 0;

double f(double ksi)
{
    return tan(ALPHA * sqrt(2 * (1 - ksi))) * sqrt(1 / ksi - 1) - 1;
}

double df(double ksi)
{
    return -(A * U0 * sqrt(1 / ksi - 1)) / (sqrt(2 * U0 * (1 - ksi)) * pow(cos(A * sqrt(2 * U0 * (1 - ksi))), 2)) - tan(A * sqrt(2 * U0 * (1 - ksi))) / (2 * ksi * ksi * sqrt(1 / ksi - 1));
}

double dichotomy(double a, double b)
{
    if (abs(b - a) > epsilon)
    {
        if (f(a) * f((a + b) / 2) <= 0)
        {
            b = (a + b) / 2;
        }
        else
        {
            a = (a + b) / 2;
        }
        // cout << a << " " << b << endl;
        n1++;
        dichotomy(a, b);
    }
    else
    {
        // cout << "a = " << a << endl;
        return (a + b) / 2;
    }
}

double simpleItterations(double x)
{
    double lambda = -0.25;
    double x_next = x;

    do
    {
        x = x_next;
        x_next = x - lambda * f(x);
        n2++;
    } while (abs(x_next - x) > epsilon);

    return (x_next + x) / 2;
}

double newton(double x)
{
    double lambda;
    double x_next = x;

    do
    {
        x = x_next;
        lambda = 1 / df(x);
        x_next = x - lambda * f(x);
        n3++;
    } while (abs(x_next - x) > epsilon);

    return (x_next + x) / 2;
}

int main()
{

    cout << "result (dichotomy) = " << dichotomy(E1, 1) << endl;
    cout << "n (dichotomy) = " << n1 << "\n"
         << endl;

    cout << "result (simpleItter) = " << simpleItterations(0.5) << endl;
    cout << "n (simpleItter) = " << n2 << "\n"
         << endl;

    cout << "result (Newton) = " << newton(0.5) << endl;
    cout << "n (Newton) = " << n3 << "\n"
         << endl;

    return 0;
}