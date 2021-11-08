import numpy as np
import matplotlib.pyplot as plt


def polynom(x, y, n, t):
    P_n = 0
    for i in range(n + 1):
        li_x = 1
        li_xi = 1
        for j in range(n + 1):
            if i != j:
                li_x *= t - x[j]
                li_xi *= x[i] - x[j]
        P_n += y[i] * li_x / li_xi
    return P_n


def main():
    n = 12
    x_i = [(1 + k / n) for k in range(n + 1)]
    y_i = [np.log(x) for x in x_i]
    s = 10 ** 2
    x = [0.5 + i / s * 20 for i in range(s + 1)]
    lagrange = [polynom(x_i, y_i, n, t) for t in x]
    log = [np.log(t) for t in x]
    result = list(map(lambda x,y: x-y, lagrange, log))
    plt.plot(x, result, 'r')
    # plt.plot(x_i, y_i, 'bo', markersize=4)
    plt.ylim(-100 / 4, 1)
    plt.xlim(1 / 2, 40 / 2)
    plt.show()


if __name__ == '__main__':
    main()