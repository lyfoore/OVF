import numpy as np
import matplotlib.pyplot as plt

def f(t,x,y):
    a = 998
    b = 1998
    c = -999
    d = -1999
    return [a*x+b*y, c*x+d*y]

def solution(t,x0,y0):
    t0 = t[0]
    a = (x0 + y0)/np.exp(-t0)
    b = -(x0 + 2*y0)/np.exp(-1000*t0)
    return [2*a*np.exp(-t) + b*np.exp(-1000*t), - a*np.exp(-t) - b*np.exp(-1000*t)]

def forward_Euler(t0,x0,y0,tend,n):
    h = (tend-t0)/(n-1)
    t = np.linspace(t0,tend,n)
    x = np.zeros(len(t))
    x[0] = x0
    y = np.zeros(len(t))
    y[0] = y0
    for i in range(1,n):
        x[i] = x[i-1] + h*f(t[i-1],x[i-1],y[i-1])[0]
        y[i] = y[i-1] + h*f(t[i-1],x[i-1],y[i-1])[1]
    return [t,x,y]

def reverse_Euler(t0,x0,y0,tend,n):
    a = 998
    b = 1998
    c = -999
    d = -1999
    h = (tend-t0)/(n-1)
    t = np.linspace(t0,tend,n)
    x = np.zeros(len(t))
    x[0] = x0
    y = np.zeros(len(t))
    y[0] = y0
    D = (1-h*a)*(1-h*d) - h*h*b*c
    for i in range(1,n):
        x[i] = ((1-h*d)*x[i-1] + h*b*y[i-1])/D
        y[i] = (h*c*x[i-1] + (1-h*a)*y[i-1])/D
    return [t,x,y]

#для прямого h < 0.002
x0 = 1
y0 = 1
n = 10000
t,x,y = forward_Euler(0,x0,y0,6,n)
x_sol, y_sol = solution(t,x0,y0)
t_rev,x_rev,y_rev = reverse_Euler(0,x0,y0,6,n)
plt.plot(t, -y, label = 'Forward Euler')
plt.plot(t, -y_rev, label = 'Reverse Euler')
plt.plot(t, -y_sol, label = 'Analitic solution')
plt.legend()
plt.show()