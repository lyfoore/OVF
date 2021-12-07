import math as m
import numpy as np
import matplotlib.pyplot as plt

a0=1
a1=0.002
w0 = 5.1
w1 = 25.5
N = 100

t1 = 0
t2 = 2*m.pi

def f(x):
    return a0*np.sin(w0*x) + a1*np.cos(w1*x)

def h(k, n):
    return 0.5*(1-m.cos(2*m.pi*k/n))


def DFT(f, N, t1, t2):
    x = np.linspace(t1, t2, N, endpoint = False)
    ftr = []
    w = []
    for i in range(0, N):
        ftri = 0
        for k in range(0, N) :
            ftri += f(*[x[k]])*np.exp(2*m.pi*(1j)*i*k/N)
        ftr.append(abs(ftri))
        w.append(2*m.pi*i/(t2-t1))
    return [w, ftr]

def DFT1(f, N, t1, t2):
    x = np.linspace(t1, t2, N, endpoint = False)
    ftr1 = []
    w = []
    for i in range(0, N):
        ftri1 = 0
        for k in range(0, N) :
            ftri1 += f(*[x[k]])*np.exp(2*m.pi*(1j)*i*k/N)*h(k,N)
        ftr1.append(abs(ftri1))
        w.append(2*m.pi*i/(t2-t1) - m.pi*N/(t2-t1))
    return [w, ftr1]

x = np.linspace(t1, t2, N)
[w, ftr] = DFT(f, N, t1, t2)
[w, ftr1] = DFT1(f, N, t1, t2)

plt.plot(w, ftr,'r-', label = 'Прямоугольное окно')
plt.yscale('log')
#plt.xlim(0,40)
plt.legend()
#plt.show()

plt.plot(w, ftr1, 'g-',label = 'Окно Ханна')
plt.yscale('log')
#plt.ylim(0.0001)
#plt.xlim(0,40)
plt.legend()
plt.show()


plt.plot(x, f(x), 'b-', label = 'f(t)')
plt.show()


