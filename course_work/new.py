import math as m
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

T_max = 1 * 1e5
R_min = 4
R_max = 15
D = 1e-1
T_N = 250
R_N = 33

dT = T_max / (T_N-1)
dR = (R_max-R_min) / (R_N-1)

alpha = dT * D / dR / dR 
beta = D * dT / dR / 2

T = np.linspace(0, T_max, T_N)
R = np.linspace(R_min, R_max, R_N)

U = np.full(R_N, 10.)

error = 1e-10
out = -1e2
inner = 1e2 
fix = 1e2

def gaussSeidel_tStep():
    Uold = U.copy()
    converge = False

    # outer neumann cond: (U[R_N-1 + 1][f] = U[R_N - 2][f] + 2*dR*out) where 'out' is magnitude of derivative
    U[R_N-1] = fix #(alpha*(U[R_N - 2] + 2*dR*out+U[R_N-2]) + beta/R[R_N-1]*(U[R_N - 2] + 2*dR*out - U[R_N-2]) + Uold[R_N-1])/(1+2*alpha)

    # inner neumann cond: (U[-1][f] = U[1][f] + 2*dR*inner) where 'inner' is magnitude of derivative
    U[0] = (2*alpha*(U[1] + dR*inner) - beta/R[0]*(2*dR*inner) + Uold[0])/(1+2*alpha)

    while not converge:
        U0 = U.copy()


        for r in reversed(range(1, R_N - 1)):
            # U[r] = (alpha * (U0[r+1] + U[r-1]) + beta/R[r]*(U0[r+1] - U[r-1]) + Uold[r])/(1 + 2*alpha)
            U[r] = (alpha*(U0[r+1]+U[r-1]) + beta/R[r]*(U0[r+1] - U[r-1]) + Uold[r])/(1+2*alpha)

        # print(np.max(abs(U - U0)))
        if (np.max(abs(U - U0)) < error):
            converge = True


def init():
    for indx, val in enumerate(U[int(R_N/2)]):
        U[int(R_N/2)][indx] = 10.

U_time = []

def draw_graph():
    U_analyt = fix + inner*R_min*np.log(R_max/R)
    plt.plot(R, U, label='calc')
    plt.plot(R, U_analyt, label='analyt')
    plt.legend()
    plt.show()

    plt.plot(T, U_time)
    plt.show()

def stat_time():
    eps = 1e-3
    stat_time = "unstationary"
    for i in range(1, T_N):
        if (abs(U_time[i] - U_time[i-1]) < eps):
            stat_time = T[i]
            break
    print("stationary time = ")
    print(stat_time)

def main():
    # draw()
    for t in range(0, T_N):
        gaussSeidel_tStep()
        print('STEP')
        U_time.append(U[int(R_N/2) - 1])
    stat_time()
    draw_graph()

main()
