import numpy as np
import matplotlib.pyplot as plt

def func(x):
    return np.sin(x)

def progonka(A,B,C,F): #b - главная, c - верхняя, a - нижняя, длина одинаковая но c[-1] и a[0] не используются
    l = len(A)
    for i in range(1, l):
        k = -A[i]/B[i-1]
        B[i] += k*C[i-1]
        F[i] += k*F[i-1]
        A[i] = 0
    y = np.zeros(l)
    y[-1] = F[-1]/B[-1]
    for i in range(l-2, -1, -1):  
        y[i] = (F[i]-C[i]*y[i+1])/B[i]
    return y


def analit_solution(x1,x2,a1,b1,c1,a2,b2,c2,x):
    if (a1 == 0) and (a2 == 0):
        print("Некорректные параметры")
    elif (a1 == 0):
        C1 = c1/b1 + np.cos(x1)
        C2 = (c2+a2*np.sin(x2)+b2*np.cos(x2) - C1*(a2*x2+b2))/a2
    else:
        C1 = (a2/a1*(c1+a1*np.sin(x1)+b1*np.cos(x1)) - (c2+a2*np.sin(x2)+b2*np.cos(x2)))/(a2/a1*(a1*x1+b1) - (a2*x2 + b2))
        C2 = (c1+a1*np.sin(x1)+b1*np.cos(x1) - C1*(a1*x1+b1))/a1
    return (-np.sin(x) + C1*x + C2)


def detpar(xmin,xmax,a1,b1,c1,a2,b2,c2,n):
    h = (xmax-xmin)/(n-1)
    X = np.linspace(xmin,xmax,n)[1:-1]
    ones_arr = np.ones(len(X))
    A = 1*ones_arr
    B = -2*ones_arr
    C = 1*ones_arr
    F = func(X)*h*h
    B[0] += -b1*A[0]/(a1*h-b1)
    F[0] += -c1*A[0]*h/(a1*h-b1)
    B[-1] += b2*C[-1]/(a2*h+b2)
    F[-1] += -c2*C[-1]*h/(a2*h+b2)
    return [A,B,C,X,F]



xmin = 0
xmax = np.pi
n = 1001
a1, b1, c1 = [0, 1, 0]
a2, b2, c2 = [1, 1, 1]

a,b,c,x,f = detpar(xmin,xmax,a1,b1,c1,a2,b2,c2,n)
y = progonka(a,b,c,f)

plt.plot(x,y, label = 'Progonka')
plt.plot(x,analit_solution(xmin,xmax,a1,b1,c1,a2,b2,c2,x), label = 'Analitic')
plt.xlabel('x')
plt.ylabel('y')
plt.legend()
plt.show()
