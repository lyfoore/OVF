from math import *
import numpy as np
import matplotlib.pyplot as plt

def U(x):
    return x**2/2
    
n = 1000
N = 100
x1 = -10
x2 = 10
m = 3

y0 = np.ones(n) + np.linspace(0, 1, n)
x = np.linspace(x1, x2, n, True)

h = x[1] - x[0]

a = [0]
for i in range(1, n):
    a.append(-1./(2*h**2))

b = []
for i in range(0, n):
    b.append(1./h**2+U(x[i]))
    
c = []
for i in range(0, n-1):
    c.append(-1./(2*h**2))
c.append(0.)

def progonka(a, b, c, d, n):
    y = []
    for i in range(0, n):
        y.append(0)
    for i in range(1, n):
        xi = a[i]/b[i-1]
        a[i] = 0
        b[i] -= xi * c[i-1]
        d[i] -= xi * d[i-1]
    y[n-1] = d[n-1]/b[n-1]     
    for i in range(n-2, -1, -1):
        y[i] = 1/b[i] * (d[i] - c[i]*y[i+1])
    return y

def obrIterr(Y0, A, B, C, n, N, m):
    psi = []
    toler = 1e-1
    E = []
    for j in range(0, m):
        Y = Y0.copy()
        for k in range(0, j):
            Y = Y - psi[k]*(np.inner(Y0, psi[k]))/np.linalg.norm(psi[k])
        for i in range(0, N): 
            new_Y = Y
            prev_Y = new_Y.copy()
            a1 = A.copy()
            b1 = B.copy()
            c1 = C.copy()
            Y = progonka(a1, b1, c1, prev_Y, n)
            for k in range(0, j):
                Y = Y - psi[k]*(np.inner(Y, psi[k]))/np.linalg.norm(psi[k])
            # print(np.linalg.norm(prev_Y)/np.linalg.norm(Y))
            # print(np.linalg.norm(prev_Y))
            # print(np.linalg.norm(Y))
            if abs(np.linalg.norm(prev_Y) - np.linalg.norm(Y) < toler):
                print('break')
                print(i)
                break
        E0 = np.linalg.norm(new_Y)/np.linalg.norm(Y)  
        Y /= np.linalg.norm(Y)
        E.append(E0)
        psi.append(Y)
    return [E, psi]


[E, psi] = obrIterr(y0, a, b, c, n, N, m)

for i in range(0, m):
    print (E[i])
    plt.plot(x, psi[i], label = 'psi' + str(i))

plt.legend()
plt.show()