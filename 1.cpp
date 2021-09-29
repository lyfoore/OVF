#include <iostream>
#include <cmath>
#include <math.h>

using namespace std;

int findEpsilon()
{
    double one = 1.;
    int i = 1;
    double var = pow(2., (-1) * i);

    while (one + var != 1.) {
        var = double(pow(2., (-1) * i));
        i++;
        // cout << var << endl;
    }

    cout << i - 1 << endl;

    return 0;
}

int findMaximum() {
    double max = 0;
    int n = 0;
    while (!isinf(max)) {
        max = max + pow(2, n);
        n++;
        // cout << max << endl;
    }
    return n;
}

int main() {
    findEpsilon();
    // cout << findMaximum() << endl;
    return 0;
}