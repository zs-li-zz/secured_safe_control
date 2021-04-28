# -*- coding: utf-8 -*-
"""
Created on Mon Mar 15 11:06:07 2021

@author: siddg
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
import control as ct
from qpsolvers import cvxopt_solve_qp as solve_qp
from matplotlib import animation
from scipy import signal
from scipy import linalg

pi = np.pi;
sin = np.sin;
cos = np.cos;
sinc = np.sinc;
exp = np.exp;
log = np.log;
log10 = np.log10;


def bandlimitedNoise(N):
    noise = np.random.normal(0, 1e-2, N)
    h = signal.firwin(256, 0.4)
    return signal.lfilter(h, 1, noise)[256:] - np.mean(signal.lfilter(h, 1, noise)[256:])


# initial = np.array([-10, -10, 0, 0])

initial = np.array([0, 11.856, 0, 0])

K_gain = np.array([0.447213595500539, 3.10711058638360, 128.125101649409, 127.003953329490]) * np.array(
    [10, 7, 1.5, 1.25])

epsilon = 1e-6
gamma = 65
zeta = 40

flag = 0

fps = 60
N = 5000

g = 9.8
M = 1
m = 0.1
l = 5
bx = 0.05
bt = 0.1

L = 35
q = 5

barrSwitch = np.array([3, 5, 3, 5])

# x_barr = 10

# v_barr = 10

# t_barr = 2.5*pi/18

# w_barr = 1

x_barr = 15

v_barr = 10

t_barr = 3 * pi / 18

w_barr = 5

P = np.zeros([4, 4])

P[0, 0] = -2 / (x_barr ** 2) * barrSwitch[0]

P[1, 1] = -2 / (v_barr ** 2) * barrSwitch[1]

P[2, 2] = -2 / (t_barr ** 2) * barrSwitch[2]

P[3, 3] = -2 / (w_barr ** 2) * barrSwitch[3]

# alpha = lambda x: np.sqrt(x*np.heaviside(x, 1)) - np.sqrt(-x*np.heaviside(-x, 0))

noise = bandlimitedNoise(N + 1000)

noise = noise / np.linalg.norm(noise, np.inf)

alpha = lambda x: x


def computeJacobian(fs, X_op):
    N = len(X_op)
    J = np.zeros((N, N))
    I = np.identity(N)
    epsilon = 1e-6

    for i in range(N):
        J[i, :] = (fs(X_op + epsilon * I[i, :]) - fs(X_op - epsilon * I[i, :])) / 2 / epsilon

    return np.transpose(J)


A = lambda X: np.array([[1, 0, 0, 0],
                        [0, m + M, 0, m * l * cos(X[2])],
                        [0, 0, 1, 0],
                        [0, m * l * cos(X[2]), 0, m * l ** 2]
                        ])

B = lambda X: np.array([X[1],
                        m * l * sin(X[2]) * X[3] ** 2 - bx * X[1],
                        X[3],
                        m * g * l * sin(X[2]) - bt * l * X[3]
                        ])

# F = lambda X, U: np.linalg.solve(A(X), B(X) + np.array([0, U, 0, 0]))

fs = lambda X: np.linalg.solve(A(X), B(X))

gs = lambda X: np.linalg.solve(A(X), np.array([0, 1, 0, 0]))

A_lin = computeJacobian(fs, np.zeros(4))

B_lin = gs(np.zeros(4))

C = np.array([[1, 0, 0, 0],
              [0, 0, 1, 0]])

# L = ct.lqr(A_lin.T, C.T, np.eye(4), np.eye(2))[0].T

# L = signal.place_poles(A_lin.T, C.T, np.array([-0.8972354 +0.54904994j, -0.8972354 -0.54904994j, -1.16923758, -1.99765673])).gain_matrix.T

L = signal.place_poles(A_lin.T, C.T, np.array([-1.25, -0.85, -0.85, -1.25])).gain_matrix.T

fs_lin = lambda X: np.matmul(A_lin, X)

gs_lin = lambda X: B_lin

F = lambda X, U: fs(X) + gs(X) * U

K = lambda X: np.dot(K_gain, X)

h = lambda X: (1 - (X[0] / x_barr) ** 2) * barrSwitch[0] + (1 - (X[1] / v_barr) ** 2) * barrSwitch[1] + (
        1 - (X[2] / t_barr) ** 2) * barrSwitch[2] + (1 - (X[3] / w_barr) ** 2) * barrSwitch[3]

h_grad = lambda X: np.matmul(P, X)

Lf_h = lambda X: np.matmul(h_grad(X).T, fs(X))

Lg_h = lambda X: np.matmul(h_grad(X).T, gs(X))

f_out = lambda X: np.matmul(C, X)

f_est = lambda X_hat, Y, U: np.matmul(A_dis, X_hat) + B_dis * U + np.matmul(L, (Y - np.matmul(C, X_hat)))


def CBFControl(X):
    K_desired = np.array([K(X)])

    qp_P = np.array([1 + epsilon])

    Lf = np.array([Lf_h(X)])
    Lg = np.array([Lg_h(X)])
    # hly = alpha(h(X))
    hly = np.array([h(X)])

    qp_q = -K_desired
    qp_h = np.array([Lf + gamma * alpha(hly) - zeta * np.matmul(Lg, Lg.T)]) - epsilon
    qp_G = -Lg - epsilon
    print("qp_P, qp_q, qp_G, qp_h", qp_P, qp_q, qp_G, qp_h)
    result = solve_qp(qp_P, qp_q, qp_G, qp_h)
    print("result: ", result)
    return result


def CBFControlNew(X):
    K_desired = np.array([K(X)])

    qp_P = np.array([1 + epsilon])

    Lf = np.array([Lf_h(X)])
    Lg = np.array([Lg_h(X)])
    # hly = alpha(h(X))
    hly = np.array([h(X)])

    qp_q = -K_desired
    qp_h = np.array([Lf + gamma * alpha(hly) - zeta * np.matmul(Lg, Lg.T)]) - epsilon
    qp_G = -Lg - epsilon

    return solve_qp(qp_P, qp_q, qp_G, qp_h)


def func(X, t, q):
    X_dot = np.zeros(8)
    i = int(t / Ts)
    U = CBFControl(X[:4])  # ORIGONAL 4:, use real states now
    X_dot[:4] = F(X[:4], U)
    X_dot[4:] = f_est(X[4:], f_out(X[:4]) + 0.01 * q * noise[i] * np.array([1, 1]), U)
    # U = K(X + 30*noise[i]*np.array([1, 1, 1, 1]))
    return X_dot


Q = 50
N = 5000
Ts = 1 / 500
t_list = np.arange(N) * Ts

e = np.zeros([N, 4])

U = np.zeros([N])

A_dis = linalg.expm(Ts*A_lin)
int_expAT = np.zeros([4,4])
for i in range(4):
    for j in range(4):
        expAT = lambda X, t: linalg.expm(t * A_lin)[i, j]
        int_expAT[i,j] = integrate.odeint(expAT, 0, [0, Ts], rtol=1e-20)[1]
B_dis = np.matmul(int_expAT, B_lin)


X = np.zeros([N, 4])
Y = np.zeros([N, 2])
X_est = np.zeros([N, 4])
X[0, :] = initial
Y[0, :] = f_out(X[0, :]) + 0.01 * q * noise[0] * np.array([1, 1]) # add output noise
U[0] = CBFControl(X[0, :])
#X_est at time 0 is zero
for i in range(1,N):
      # calculate control with states at last time step
    X[i, :] = np.matmul(A_dis, X[i-1, :]) + B_dis * U[i-1]
    Y[i, :] = f_out(X[i, :]) + Ts * 0.01 * q * noise[i] * np.array([1, 1]) # add output noise
    X_est[i, :] = f_est(X[i - 1, :], Y[i, :], U[i-1])
    U[i] = CBFControl(X[i, :])
    print("step num: ", i)
    print("X_RESULT: ", X[i, :])

e = X - X_est
## ORIGINAL CODE
# X = np.zeros([N, 4])
# X[0, :] = initial
# for i in range(N-1):
#     U[i] = CBFControl(X[i, :]) # calculate control with states at last time step
#     X[i+1, :] = np.matmul(A_dis, X[i, :]) + B_dis*U[i]
#     print("X_RESULT: ", X[i, :])


# for q in range(Q):
#     X[:, :, q] = integrate.odeint(func, initial, t, args = (q, ))

# u_norm = np.zeros(Q)
# e_norm = np.zeros([Q, 4])
# n_norm = np.arange(Q)*0.01

# for q in range(Q):
# for i in range(N):
#     U[i] = CBFControl(X[i, 4:])

#     u_norm[q] = np.linalg.norm(U[:, q], np.inf)


# for q in range(Q):
#     e[:, :, q] = X[:, 4:, q] - X[:, :4, q]
#     for j in range(4):
#         e_norm[q, j] = np.linalg.norm(e[:, j, q], np.inf)


plt.figure()

plt.subplot(211)
plt.plot(t_list, X)
plt.xlabel('time in seconds')
plt.title('Actual system states')
plt.grid()
# plt.subplot(212)
# plt.plot(t_list, U)
# plt.xlabel('time in seconds')
# plt.title('Control Input')
# plt.grid()

plt.subplot(212)

plt.plot(t_list, e)
plt.ylim([-0.1, 0.1])
plt.xlabel('time in seconds')

plt.title('Estimated ERROR')

plt.grid()

plt.savefig("states.png")

# plt.subplot(121)

# plt.plot(n_norm, e_norm[:, 0], label = '$x - \hat{x}||$')
# plt.plot(n_norm, e_norm[:, 1], label = '$\dot{x} - \dot{\hat{x}}$')
# plt.plot(n_norm, e_norm[:, 2], label = '$\Theta - \hat{\Theta}$')
# plt.plot(n_norm, e_norm[:, 3], label = '$\dot{\Theta} - \dot{\hat{\Theta}}$')

# plt.grid()

# plt.xlabel('$||n||$')

# plt.ylabel('$||e||$')

# plt.title('Infinity Norm of State Estimator Error')

# plt.legend()

# plt.subplot(122)

# plt.plot(n_norm, u_norm)

# plt.grid()

# plt.xlabel('$||n||$')

# plt.ylabel('$||u||$')

# plt.title('Infinity Norm of Controller Effort')

h_test = (1 - (X[:, 0] / x_barr) ** 2) * barrSwitch[0] + (1 - (X[:, 1] / v_barr) ** 2) * barrSwitch[1] + (
        1 - (X[:, 2] / t_barr) ** 2) * barrSwitch[2] + (1 - (X[:, 3] / w_barr) ** 2) * barrSwitch[3]

# xc = X[:, 0]
# theta = X[:, 2]

# xp = xc + l*sin(theta)
# yp = l*cos(theta)

plt.figure()
plt.plot(t_list, h_test)
plt.axhline(y=0, c='r')
plt.grid()
plt.xlabel('Time(s)', fontsize=18)
plt.ylabel('Barrier Zeroing Function $h(X)$', fontsize=18)
plt.savefig("test_CBF.png")

# Writer = animation.writers['ffmpeg']
# writer = Writer(fps=60, metadata=dict(artist='Me'), bitrate=2400)

# fig = plt.figure(figsize=(18, 3.5), dpi=80, )
# ax = plt.axes(xlim=(-L, L), ylim=(-1.1*l, 1.1*l))
# line1, = ax.plot([], [], lw=2, marker = 'o', c = 'b')
# line2, = ax.plot([], [], marker = 's', linestyle = 'None', c = '#02b03c', markersize = 10)
# plt.gca().set_aspect('equal', adjustable='box'); plt.tight_layout(); plt.axhline(y = 0, c = 'k', linewidth = 2)
# plt.grid()
# plt.title('Pendulum on a cart with LQR and linearized CBF')

# def init():
#     line1.set_data([], [])
#     line2.set_data([], [])
#     return line1, line2,

# # animation function.  This is called sequentially
# def animate(i):
#     line1.set_data([xc[i], xp[i]], [0, yp[i]])
#     line2.set_data([xc[i]], [0])
#     plt.plot([-L, L], [0, 0], lw = 1.5, c = 'k')
#     return line1, line2,

# # call the animator.  blit=True means only re-draw the parts that have changed.
# anim = animation.FuncAnimation(fig, animate, init_func=init,
#                                 frames=N, interval=1000/fps, blit=True)
# plt.tight_layout()
# #anim.save('pendNoise1.mp4', writer=writer)