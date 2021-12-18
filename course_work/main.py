import math as m
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

T_max = 1* 1e2
R_min = 0.4
R_max = 1
D = 1e-3
T_N = 250
R_N = 21
F_N = 10

dT = T_max / T_N
dR = R_max / R_N
dF = 2 * np.pi / F_N

alpha = dT * D / dR / dR
beta = D * dT / dR
gamma = D * dT / dF / dF * 0

T = np.linspace(0, T_max, T_N)
R = np.linspace(R_min, R_max, R_N)
F = np.linspace(0, 2 * np.pi, F_N)

U = np.full((R_N, F_N), 0.)

error = 1e-6
out = 1e4   # + for out
inner = -1e4  # - for in
def gaussSeidel_tStep():
    Uold = U.copy()
    converge = False
    while not converge:
        U0 = U.copy()

        # outer neumann cond: (U[R_N-1 + 1][f] = U[R_N - 2][f] - 2*dR*out) where 'out' is magnitude of derivative
        # U[R_N-1][0] = (alpha*(U[R_N - 2][0] - 2*dR*out + Uold[R_N-2][0]) + beta/R[R_N-1]*(U[R_N - 2][0] - 2*dR*out) + gamma/R[R_N-1]/R[R_N-1]*(U[R_N-1][1] + Uold[R_N-1][F_N-1]) + Uold[R_N-1][0])/(1+2*alpha+beta/R[R_N-1]+2*gamma/R[R_N-1]/R[R_N-1])
        for f in range(0, F_N):
            if f == 0:
                U[R_N-1][f] = (alpha*(U[R_N - 2][f] - 2*dR*out + U[R_N-2][f]) + beta/R[R_N-1]*(U[R_N - 2][f] - 2*dR*out - U[R_N-2][f]) + gamma /R[R_N-1]/R[R_N-1]*(U[R_N-1][f+1] + U[R_N-1][f-1]) + Uold[R_N-1][f])/(1+2*alpha+2/R[R_N-1]/R[R_N-1]*gamma)
            U[R_N-1][f] = U[R_N-1][0]
        # U[R_N-1][F_N - 1] = (alpha*(U[R_N - 2][F_N - 1] - 2*dR*out + Uold[R_N-2][F_N - 1]) + beta/R[R_N-1]*(U[R_N - 2][F_N - 1] - 2*dR*out) + gamma/R[R_N-1]/R[R_N-1]*(U[R_N-1][0] + Uold[R_N-1][F_N-2]) + Uold[R_N-1][F_N - 1])/(1+2*alpha+beta/R[R_N-1]+2*gamma/R[R_N-1]/R[R_N-1])

        # inner neumann cond: (U[-1][f] = U[1][f] - 2*dR*inner) where 'inner' is magnitude of derivative
        # U[0][f] = (alpha*(U[1][f] + Uold[1][f] - 2*dR*inner) + beta/R[0]*U[1][f] + gamma/R[0]/R[0]*(U[0][f+1] + Uold[0][f-1]) + Uold[0][f])/(1+2*alpha+beta/R[0]+2*gamma/R[0]/R[0])
        for f in range(0, F_N):
            if f == 0:
                U[0][f] = (alpha*(U0[1][f] + U0[1][f] - 2*dR*inner) + beta/R[0]*(U0[1][f] - U0[1][f] - 2*dR*inner) + gamma /R[0]/R[0]*(U[0][f+1] + U[0][f-1]) + Uold[0][f])/(1+2*alpha+2/R[0]/R[0]*gamma)
            U[0][f] = U[0][0]
        # U[0][F_N - 1] = (alpha*(U[R_N - 2][F_N - 1] - 2*dR*out + Uold[R_N-2][F_N - 1]) + beta/R[R_N-1]*(U[R_N - 2][F_N - 1] - 2*dR*out) + gamma/R[R_N-1]/R[R_N-1]*(U[R_N-1][0] + Uold[R_N-1][F_N-2]) + Uold[R_N-1][F_N - 1])/(1+2*alpha+beta/R[R_N-1]+2*gamma/R[R_N-1]/R[R_N-1])


        for r in range(1, R_N - 1):
            # complete the periodicaly condition for angle f. U[r][0]=U[r][F_N-1]
            # U[r][0] = (alpha*(U[r+1][0] + Uold[r-1][0]) + beta/R[r]*U[r+1][0] + gamma/R[r]/R[r]*(U[r][1] + Uold[r][F_N - 1]) + Uold[r][0])/(1+2*alpha+beta/R[r]+2*gamma/R[r]/R[r])
            for f in range(0, F_N):
                if f == 0:
                    U[r][f] = (alpha*(U0[r+1][f] + U[r-1][f]) + beta/R[r]*(U0[r+1][f] - U[r-1][f]) + gamma /R[r]/R[r]*(U0[r][f+1] + U[r][f-1]) + Uold[r][f])/(1+2*alpha+2/R[r]/R[r]*gamma)
                U[r][f] = U[r][0]
            # U[r][F_N - 1] = (alpha*(U[r+1][F_N - 1] + Uold[r-1][F_N - 1]) + beta/R[r]*U[r+1][F_N - 1] + gamma/R[r]/R[r]*(U[r][0] + Uold[r][F_N - 2]) + Uold[r][F_N - 1])/(1+2*alpha+beta/R[r]+2*gamma/R[r]/R[r])

        # print(np.max(abs(U - U0)))
        if (np.max(abs(U - U0)) < error):
            converge = True



def init():
    for indx, val in enumerate(U[int(R_N/2)]):
        U[int(R_N/2)][indx] = 10.
        # U[int(R_N/2)+1][indx] = 10.


def draw():
    Tmax = np.max(U)
    Tmin = -Tmax
    fig = plt.figure()
    ax = fig.add_subplot(122, polar=True)
    bx = fig.add_subplot(121, polar=True)
    ax.set_yticklabels([])
    bx.set_yticklabels([])

    plt.setp([ax, bx], rorigin=0, rmin=0.0, rmax=R_max)

    ctf = ax.contourf(F, R, U, 50, cmap=cm.seismic, vmin=Tmin, vmax=Tmax)
    ctf_old = bx.contourf(F, R, U_init, 500, cmap=cm.seismic, vmin=Tmin, vmax=Tmax)
    cbar_ax = fig.add_axes([0.9, 0.15, 0.03, 0.7])
    fig.subplots_adjust(right=0.85, left=0.08)
    fig.colorbar(ctf, cax=cbar_ax)
    plt.show()

def main():
    # draw()
    for t in range(0, T_N):
        gaussSeidel_tStep()
        print('STEP')
    print(U[int(R_N/2)][int(F_N/2)])
    print(U[int(R_N/2)-1][int(F_N/2)])
    print(U[int(R_N/2)+1][int(F_N/2)])
    draw()


# init()
U_init = U.copy()
main()
