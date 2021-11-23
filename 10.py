import numpy as np
import matplotlib.pyplot as plt

def left_border(t):
    return 0

def right_border(t):
    return 0

def start_T(x, L):
    return x*(1 - x/L)**2

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

def par(prev_y, t, h, dt, Amp):
    l = len(prev_y)
    ones_arr = np.ones(l)
    A = 1*ones_arr
    B = -2*(1+h*h/dt/Amp)*ones_arr
    C = 1*ones_arr
    F = np.zeros(l)
    
    F[0] = -2*(h*h/dt/Amp-1)*prev_y[0] - prev_y[1] - left_border(t)
    for i in range(1,l-1):
        F[i] = -2*(h*h/dt/Amp-1)*prev_y[i] - prev_y[i+1] - prev_y[i-1]
    F[-1] = -2*(h*h/dt/Amp-1)*prev_y[-1] - prev_y[-2] - right_border(t)
    
    return [A,B,C,F]

def derivative(x, y):
    l = len(y)-1
    diff_y = np.zeros(l)
    new_x = np.zeros(l)
    for i in range(l):
        diff_y[i] = (y[i+1]-y[i])/(x[i+1]-x[i])
        new_x[i] = (x[i+1]+x[i])/2
    return [new_x, diff_y]

xmin = 0
xmax = 1
n = 101
T = 0.5
nT = 1001
Amp = 1

arr_t = np.linspace(0,T,nT)
arr_Tmax = np.linspace(0,T,nT)
dt = T/(nT-1)
L = xmax - xmin
h = L/(n-1)
X = np.linspace(xmin, xmax, n)
Y = start_T(X, L)
for i in range(len(arr_t)):
    t = arr_t[i]
    A,B,C,F = par(Y, t, h, dt, Amp)
    Y = progonka(A,B,C,F)
    if (1==0):
        plt.xlabel('x')
        plt.ylabel("T")
        plt.plot(X, Y)
        plt.show()
    arr_Tmax[i] = max(Y)

plt.plot(X, Y)
plt.xlabel('X')
plt.ylabel("T")
plt.show()
print (Y[0], Y[-1])

X1, Y1 = derivative(X, Y)
plt.plot(X1, Y1)
plt.xlabel('X')
plt.ylabel("T'")
plt.show()
    
plt.plot(np.log(arr_t[1:-1]), np.log(arr_Tmax[1:-1]))
plt.xlabel('ln(t)')
plt.ylabel("ln(Tmax)")
plt.show()


