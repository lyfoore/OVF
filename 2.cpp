#include <iostream>
#include <cmath>

#define A 1
#define U0 1
#define epsilon 2e-16

using namespace std;

int n = 0;

double f(double ksi) {
    return tan(A*sqrt(2*U0*(1 - ksi)))*sqrt(1/ksi - 1) - 1;
}

double dichotomy(double a, double b) {
    if (abs(b - a) > epsilon) {
        if (f(a) * f((a + b) / 2) <= 0) {
            b = (a + b) / 2;
        }
        else {
            a = (a + b) / 2;
        }
        // cout << a << " " << b << endl;
        n++;
        dichotomy(a, b);
    }
    else {
        // cout << "a = " << a << endl;
        return a;
    }
}

double simpleItterations (double x_0) {
    double x;
    double x_next;
    while (abs(x - x_next) > epsilon) {
        x = x_next;
        x_next = f(x_0) + x_0;
        x_0 = x_next;
        cout << x_next << endl;
    }
    return x_next;
}

int main() {
    // double result = dichotomy(0, 1);

    // cout << "result (dichotomy) = " << dichotomy(0, 1) << endl;
    // cout << "n = "<< n << endl;

    cout << simpleItterations(0.5);
    
    return 0;
}