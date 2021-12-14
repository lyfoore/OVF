import math as m
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

T_max = 1
R_min = 0.5
R_max = 1
D = 1e1
T_N = 10
R_N = 100
F_N = 100

dT = T_max / T_N
dR = R_max / R_N
dF = 2 * np.pi / F_N

alpha = dT * D / dR / dR
beta = D * dT / dR
gamma = D * dT / dF / dF

T = np.linspace(0, T_max, T_N)
R = np.linspace(R_min, R_max, R_N)
R_full = np.linspace(0, R_max, R_N)
F = np.linspace(0, 2 * np.pi, F_N)

U = np.full((R_N, F_N), 0.)

error = 1e-5
out = -1e2
def gaussSeidel_tStep():
    Uold = U.copy()
    converge = False
    while not converge:
        U0 = U.copy()

        # neumann cond: (U[R_N-1 + 1][f] = U[R_N - 2][f] - 2*dR*out) where out is magnitude of derivative
        for f in range(1, F_N - 1):
            U[R_N-1][f] = (alpha*(U[R_N - 2][f] - 2*dR*out + Uold[R_N-2][f]) + beta/R[R_N-1]*(U[R_N - 2][f] - 2*dR*out) + gamma/R[R_N-1]/R[R_N-1]*(U[R_N-1][f+1] + Uold[R_N-1][f-1]) + Uold[R_N-1][f])/(1+2*alpha+beta/R[R_N-1]+2*gamma/R[R_N-1]/R[R_N-1])

        for r in range(1, R_N - 1):
            # complete the periodicaly condition for angle f. U[r][0]=U[r][F_N-1]
            U[r][0] = (alpha*(U[r+1][0] + Uold[r-1][0]) + beta/R[r]*U[r+1][0] + gamma/R[r]/R[r]*(U[r][1] + Uold[r][R_N - 1]) + Uold[r][0])/(1+2*alpha+beta/R[r]+2*gamma/R[r]/R[r])
            for f in range(1, F_N - 1):
                U[r][f] = (alpha*(U[r+1][f] + Uold[r-1][f]) + beta/R[r]*U[r+1][f] + gamma/R[r]/R[r]*(U[r][f+1] + Uold[r][f-1]) + Uold[r][f])/(1+2*alpha+beta/R[r]+2*gamma/R[r]/R[r])
            U[r][R_N - 1] = (alpha*(U[r+1][R_N - 1] + Uold[r-1][R_N - 1]) + beta/R[r]*U[r+1][R_N - 1] + gamma/R[r]/R[r]*(U[r][0] + Uold[r][R_N - 2]) + Uold[r][R_N - 1])/(1+2*alpha+beta/R[r]+2*gamma/R[r]/R[r])

        # print(np.max(abs(U - U0)))
        if (np.max(abs(U - U0)) < error):
            converge = True



def init():
    for indx, val in enumerate(U[int(R_N/2)]):
        U[int(R_N/2)][indx] = 10.
        U[int(R_N/2)+1][indx] = 10.


def draw():
    Tmax = np.max(U)
    Tmin = 0.
    fig = plt.figure()
    ax = fig.add_subplot(122, polar=True)
    bx = fig.add_subplot(121, polar=True)
    ax.set_yticklabels([])
    bx.set_yticklabels([])
    # ax.set_rscale('log')
    # ax.set(rmax=1.0, rmin=0.0)
    plt.setp([ax, bx], rorigin=0, rmin=5, rmax=10)

    ctf = ax.contourf(F, R, U, 50, cmap=cm.hot, vmin=Tmin, vmax=Tmax)
    ctf_old = bx.contourf(F, R, U_init, 50, cmap=cm.hot, vmin=Tmin, vmax=Tmax)
    cbar_ax = fig.add_axes([0.9, 0.15, 0.03, 0.7])
    fig.subplots_adjust(right=0.85, left=0.08)
    fig.colorbar(ctf_old, cax=cbar_ax)
    plt.show()

def main():
    # draw()
    for t in range(0, T_N):
        gaussSeidel_tStep()
        print('STEP')
    draw()


init()
U_init = U.copy()
main()
