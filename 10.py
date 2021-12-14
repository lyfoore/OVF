import numpy as np
import matplotlib.pyplot as plt

def left_border(t):
    return 0

def right_border(t):
    return 0

def start_T(x, L):
    return x*(1 - x/L)**2

def progonka(A,B,C,F): 
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
    
analyt_Tmax = np.linspace(0,T,nT)
analyt_Tmax = 4 / (np.pi ** 3 ) * np.exp(-np.pi**2 * analyt_Tmax)

plt.plot((arr_t[1:-1]), (analyt_Tmax[1:-1]))
plt.plot((arr_t[1:-1]), (arr_Tmax[1:-1]))
plt.xlabel('t')
plt.ylabel("Tmax")
plt.show()


