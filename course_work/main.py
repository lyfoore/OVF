import math as m
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

T_max = 1
R_max = 1
D = 10.
T_N = 100
R_N = 100
F_N = 100

dT = T_max / T_N
dR = R_max / R_N
dF = 2 * np.pi / F_N

alpha = dT / D
gamma = dT / D / dR / dR
beta = 1 / D / dT

T = np.linspace(0, T_max, T_N)
R = np.linspace(0, R_max, R_N)
F = np.linspace(0, 2 * np.pi, F_N)

U = np.full((T_N, R_N, F_N), 0.)
error = 1e-5

def gaussSeidel_tStep(T):   # where T is index
    converge = False
    while not converge:
        U0 = U.copy()
        for r in range(1, R_N - 1):
            for f in range(1, F_N - 1):
                    U[T][r][f] = gamma * (U0[T][r+1][f] + U0[T][r-1][f]) + alpha / R[r] / dR * U[T][r+1][f] + beta * (U0[T][r][f+1] + U0[T][r][f-1]) / R[r] / R[r] + U0[T-1][r][f]
        print(max(max(abs(U - U0))))
        if (max(max(abs(U - U0))) < error):
            converge = True



def init():
    for indx, val in enumerate(U[0][int(R_N/2)]):
        U[0][int(R_N/2)][indx] = 1.


def draw(t):
    ax = plt.subplot(111, polar=True)
    ax.set_yticklabels([])
    ctf = ax.contourf(F, R, U[t], cmap=cm.hot)
    plt.colorbar(ctf)
    plt.show()

def main():
    # for t in range(0, T_N + 1):
    gaussSeidel_tStep(1)
    print("DONE")


init()
main()

# gaussSeidel()
# draw(0)
# draw(1)
# print(U[0])
# for x in range(10 + 1):
#     print(1/10*x)


# function out = gs_trans_implicit(dt,dx,dy,nt,nx,ny,T,Told,implicit,x,y)

# %terms used in the equation
# gs_iter=0;
# alpa=0.5;
# k1=alpa*(dt/dx^2);
# k2=alpa*(dt/dy^2);
# term=1/(1+2*k1+2*k2);
# error=1e-4;
# tol=1e-5;

# for k = 1:length(nt)#time loop
#     while(error>tol)#convergence_loop
#         for i = 2:nx-1
#             for j = 2:ny-1
# %                   T(i,j)=term*(Told(i,j) )+ k1*(Told(i-1,j)+Told(i+1,j)) + k2*(Told(i,j-1)+Told(i,j+1));
#                 T(i,j)=term*(Told(i,j) + k1*(T(i-1,j)+T(i+1,j)) + k2*(T(i,j-1)+T(i,j+1)));
#     error=max(max(abs(T-Told)));
#     Told=T;#updating old values
#     gs_iter=gs_iter+1;
#     figure(implicit)
#     [C,d]=contourf(x,y,T);
#     colormap(jet)
#     colorbar
#     xlabel('X-axis')
#     ylabel('Y-axis')
#     clabel(C,d)
#     heading=['               Iteration number = ', num2str(gs_iter);,'when Gauss Seidel technique is applied' ];
#     title(heading)
#     pause(0.005)
#     endwhile
# iteration_value=fprintf('No. of iterations =%d when Gauss Seidel method is applied implicitlyn',gs_iter);
# out=figure(implicit);

# end 