import math as m
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

T_max = 1
R_max = 1
D = 0.1
T_N = 50
R_N = 100
F_N = 100

dT = T_max / T_N
dR = R_max / R_N
dF = 2 * np.pi / F_N

alpha = dT * D / dR / dR
beta = D * dT / dR
gamma = D * dT / dF / dF

T = np.linspace(0, T_max, T_N)
R = np.linspace(0, R_max, R_N)
F = np.linspace(0, 2 * np.pi, F_N)

U = np.full((R_N, F_N), 0.)

error = 1e-5

def gaussSeidel_tStep():
    Uold = U.copy()
    converge = False
    while not converge:
        U0 = U.copy()
        for r in range(1, R_N - 1):
            U[r][0] = (alpha*(U[r+1][0] + Uold[r-1][0]) + beta/R[r]*U[r+1][0] + gamma/R[r]/R[r]*(U[r][1] + Uold[r][R_N - 1]) + Uold[r][0])/(1+2*alpha+beta/R[r]+2*gamma/R[r]/R[r])
            for f in range(1, F_N - 1):
                U[r][f] = (alpha*(U[r+1][f] + Uold[r-1][f]) + beta/R[r]*U[r+1][f] + gamma/R[r]/R[r]*(U[r][f+1] + Uold[r][f-1]) + Uold[r][f])/(1+2*alpha+beta/R[r]+2*gamma/R[r]/R[r])
            U[r][R_N - 1] = (alpha*(U[r+1][R_N - 1] + Uold[r-1][R_N - 1]) + beta/R[r]*U[r+1][R_N - 1] + gamma/R[r]/R[r]*(U[r][0] + Uold[r][R_N - 2]) + Uold[r][R_N - 1])/(1+2*alpha+beta/R[r]+2*gamma/R[r]/R[r])

        # print(np.max(abs(U - U0)))
        if (np.max(abs(U - U0)) < error):
            converge = True



def init():
    for indx, val in enumerate(U[int(R_N/2)]):
        U[int(R_N/2)][indx] = 1.
        U[int(R_N/2)+1][indx] = 1.


def draw():
    ax = plt.subplot(111, polar=True)
    ax.set_yticklabels([])
    ctf = ax.contourf(F, R, U, cmap=cm.hot)#, vmin=0., vmax=1.)
    plt.colorbar(ctf)
    plt.show()

def main():
    draw()
    for t in range(0, int((T_N + 1))):
        gaussSeidel_tStep()
        # print('STEP')
    draw()


init()
main()
