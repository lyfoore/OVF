from math import *
import numpy as np
import matplotlib.pyplot as plt

a0=1
a1=0.2
w0 = 5.1
w1 = 25.1
N = 106

t1 = 0
t2 = 2*np.pi
T = t2 - t1

def f(x):
    return a0*np.sin(w0*x) + a1*np.sin(w1*x)

def h(k, n):
    return 0.5*(1-cos(2*np.pi*k/n))


def DFT(f):
    f_1=np.empty(len(f))
    for k in range(len(f)):
        sum=0
        for i in range(len(f)):
	        sum+=f[i]*np.exp((2*np.pi*1j*i*k)/len(f))
        f_1[k]=sum/(len(f)**0.5)
    return f_1

def DFT_HANN(f):
    f_1=np.empty(len(f))
    for k in range(len(f)):
        sum=0
        for i in range(len(f)):
            w=0.5*(1-np.cos(2*np.pi*i/(len(f))))
            sum+=w*f[i]*np.exp((2*np.pi*complex(0,1)*i*k)/len(f))
        f_1[k]=sum/(len(f)**0.5)
    return f_1

t=np.linspace(0,T,N)
w=np.linspace(0,N,N)*2*np.pi/T

signal = f(t)

ftr = DFT(signal)
ftr1 = DFT_HANN(signal)

# print(w)

# fix
# for i in ftr:

# print(ftr1)
# print(len(ftr1))
# ftr2 = ftr1[50:] + ftr1[:50]
# ftr3 = ftr[50:] + ftr[:50]
# ftr2 = ftr2.append(ftr[0 : 50])


plt.plot(w, abs(ftr.real),'r-', label = 'Rectangle Window')
plt.yscale('log')
# plt.xlim(0,50)
plt.legend()
#plt.show()

plt.plot(w, abs(ftr1.real), 'g-',label = 'Hann Window')
plt.yscale('log')
#plt.ylim(0.0001)
# plt.xlim(0,N/2)
plt.legend()
plt.show()

# plt.plot(t, f(t), 'b-', label = 'f(t)')
# plt.show()
