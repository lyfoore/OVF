import numpy as np
import matplotlib.pyplot as plt

def f(t,x,y):
    a = d = 10
    b = c = 2
    return [a*x - b*x*y, c*x*y - d*y]

def solver_Runge2(t0,x0,y0,tend,n):
    h = (tend-t0)/(n-1)
    t = np.linspace(t0,tend,n)
    x = np.zeros(n)
    y = np.zeros(n)
    x[0] = x0
    y[0] = y0
    for i in range(1,n):
        x[i] = x[i-1] + h/4*(f(t[i-1],x[i-1],y[i-1])[0] + 3*f(t[i-1]+2/3*h, x[i-1]+2/3*h*f(t[i-1],x[i-1],y[i-1])[0], y[i-1]+2/3*h*f(t[i-1],x[i-1],y[i-1])[1])[0])
        y[i] = y[i-1] + h/4*(f(t[i-1],x[i-1],y[i-1])[1] + 3*f(t[i-1]+2/3*h, x[i-1]+2/3*h*f(t[i-1],x[i-1],y[i-1])[0], y[i-1]+2/3*h*f(t[i-1],x[i-1],y[i-1])[1])[1])
    return [t,x,y]

t,x,y = solver_Runge2(0,5,5.1,5.1,1500)
plt.plot(x, y, label = 'Phase curve')
plt.legend()
plt.show()
plt.plot(t, x, label = 'x')
plt.plot(t, y, label = 'y')
plt.legend()
plt.show()