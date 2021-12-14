import numpy as np
import matplotlib.pyplot as plt

def f(x,y):
    return np.sin(y)

def y_analitic(x):
    return np.exp(-x)

def solver_Euler(x0,y0,xend,n):
    h = (xend-x0)/(n-1)
    x = np.linspace(x0,xend,n)
    y = np.zeros(len(x))
    y[0] = y0
    for i in range(1,n):
        y[i] = y[i-1] + h*f(x[i-1],y[i-1])
    return [x,y]

def solver_Runge2(x0,y0,xend,n):
    h = (xend-x0)/(n-1)
    x = np.linspace(x0,xend,n)
    y = np.zeros(len(x))
    y[0] = y0
    for i in range(1,n):
        y[i] = y[i-1] + h/4*(f(x[i-1],y[i-1]) + 3*f(x[i-1]+2/3*h, y[i-1]+2/3*h*f(x[i-1],y[i-1])))
    return [x,y]

def solver_Runge4(x0,y0,xend,n):
    h = (xend-x0)/(n-1)
    x = np.linspace(x0,xend,n)
    y = np.zeros(len(x))
    y[0] = y0
    for i in range(1,n):
        k1 = f(x[i-1], y[i-1])
        k2 = f(x[i-1]+h/2, y[i-1]+h/2*k1)
        k3 = f(x[i-1]+h/2, y[i-1]+h/2*k2)
        k4 = f(x[i-1]+h, y[i-1]+h*k3)
        y[i] = y[i-1] + h/6*(k1 + 2*k2 + 2*k3 + k4)
    return [x,y]

def drawer(x0,y0,xend,n):
    x,y1 = solver_Euler(x0,y0,xend,n)
    plt.plot(x, y1, label = 'Euler')
    x,y2 = solver_Runge2(x0,y0,xend,n)
    plt.plot(x, y2, label = 'Runge2')
    x,y3 = solver_Runge4(x0,y0,xend,n)
    plt.plot(x, y3, label = 'Runge4')
    y = y_analitic(x)
    plt.plot(x,y,label = 'Analitic')
    plt.legend()
    plt.title("n = "+str(n))
    plt.show()
    
    plt.plot(x, np.abs((y1-y)/y), label = 'Euler')
    plt.plot(x, np.abs((y2-y)/y), label = 'Runge2')
    plt.plot(x, np.abs((y3-y)/y), label = 'Runge4')
    plt.legend()
    plt.title("n = "+str(n))
    plt.show()

drawer(0, 1, 3, 50)

def max_error(x0,y0,xend,n,func):
    x,y = func(x0,y0,xend,n)
    y_an = y_analitic(x)
    return(max(np.abs((y-y_an)/y_an)))

n = [5,10,20,40,80]
errors_Euler = np.zeros(len(n))
errors_Runge2 = np.zeros(len(n))
errors_Runge4 = np.zeros(len(n))
for i in range(len(n)):
    errors_Euler[i] = max_error(0,1,3,n[i],solver_Euler)
    errors_Runge2[i] = max_error(0,1,3,n[i],solver_Runge2)
    errors_Runge4[i] = max_error(0,1,3,n[i],solver_Runge4)

n_log = np.log(n)
errors_Euler_log = np.log(errors_Euler)
errors_Runge2_log = np.log(errors_Runge2)
errors_Runge4_log = np.log(errors_Runge4)
    
print("Euler n^: " + str((errors_Euler_log[-1]-errors_Euler_log[0])/(n_log[-1]-n_log[0])*(-1)))
print("Runge2 n^: " + str((errors_Runge2_log[-1]-errors_Runge2_log[0])/(n_log[-1]-n_log[0])*(-1)))
print("Runge4 n^: " + str((errors_Runge4_log[-1]-errors_Runge4_log[0])/(n_log[-1]-n_log[0])*(-1)))

# plt.plot(n_log, errors_Euler_log)
# plt.plot(n_log, errors_Runge2_log)
# plt.plot(n_log, errors_Runge4_log) 
plt.show()