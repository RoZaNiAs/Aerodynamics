# -*- coding: utf-8 -*-
"""
Created on Mon Apr  7 16:56:38 2025

@author: mdeto
"""
import numpy as np
import matplotlib.pyplot as plt

f, xf, t = 0.02, 0.3, 0.19
N = 220
ca, deltaa = 0.3, 5 * (np.pi / 180)
alpha = 10 * np.pi / 180
timon = False

beta = np.linspace(0, 2 * np.pi, N + 1)
x = 0.5 * (1 + np.cos(beta))

if timon:
    y = np.empty_like(x)
    for i in range(len(x)):
        if x[i] <= xf:
            y[i] = (f / (xf**2)) * (2 * xf * x[i] - x[i]**2)
        elif x[i] <= (1 - ca):
            y[i] = (f / (1 - xf)**2) * ((1 - 2 * xf) + 2 * xf * x[i] - x[i]**2)
        else:
            y[i] = ((f / (1 - xf)**2) * ((1 - 2 * xf) + 2 * xf * x[i] - x[i]**2)) - (x[i] - (1 - ca)) * np.tan(deltaa)
else:
    y = np.where(
        x <= xf,
        (f / (xf**2)) * (2 * xf * x - x**2),
        (f / (1 - xf)**2) * ((1 - 2 * xf) + 2 * xf * x - x**2)
    )

yt = 5 * t * (0.2969 * np.sqrt(x) - 0.126 * x - 0.3516 * x**2 + 0.2843 * x**3 - 0.1015 * x**4)
yext = y + yt
yint = y - yt
mid = N // 2 + 1
yp = np.concatenate((yint[:mid], yext[mid:]))
xp = np.concatenate((x[:mid], x[mid:]))

x_diff = xp[1:] - xp[:-1]
y_diff = yp[1:] - yp[:-1]
theta = np.arctan2(y_diff, x_diff)
cuerda = x_diff * np.cos(theta) + y_diff * np.sin(theta)

x_pc = 0.5 * (xp[:-1] + xp[1:])
y_pc = 0.5 * (yp[:-1] + yp[1:])

cos_theta = np.cos(theta)
sin_theta = np.sin(theta)

dx = x_pc[:, None] - xp[None, :-1]
dy = y_pc[:, None] - yp[None, :-1]

x_ij = dx * cos_theta + dy * sin_theta
z_ij = -dx * sin_theta + dy * cos_theta

r1 = np.sqrt(x_ij**2 + z_ij**2)
r2 = np.sqrt((x_ij - cuerda)**2 + z_ij**2)

theta1 = np.arctan2(z_ij, x_ij)
theta2 = np.arctan2(z_ij, x_ij - cuerda)

log_term = np.log(r2 / r1) / (2 * np.pi)
angle_diff = (theta2 - theta1) / (2 * np.pi)

u_q = -log_term
w_q = angle_diff
u_g = angle_diff
w_g = log_term

np.fill_diagonal(w_q, 0.5)
np.fill_diagonal(u_g, 0.5)

uq_ij = u_q * cos_theta - w_q * sin_theta
wq_ij = u_q * sin_theta + w_q * cos_theta
ug_ij = u_g * cos_theta - w_g * sin_theta
wg_ij = u_g * sin_theta + w_g * cos_theta

a_ij = np.zeros((N + 1, N + 1))
p_ij = -ug_ij * np.sin(theta[:, None]) + wg_ij * np.cos(theta[:, None])
a_ij[:N, :N] = -uq_ij * np.sin(theta[:, None]) + wq_ij * np.cos(theta[:, None])
a_ij[:N, N] = np.sum(p_ij, axis=1)
a_ij[N, :N] = (
    uq_ij[0, :] * np.cos(theta[0]) + wq_ij[0, :] * np.sin(theta[0]) +
    uq_ij[-1, :] * np.cos(theta[-1]) + wq_ij[-1, :] * np.sin(theta[-1])
)
s_ij = (
    ug_ij[0, :] * np.cos(theta[0]) + wg_ij[0, :] * np.sin(theta[0]) +
    ug_ij[-1, :] * np.cos(theta[-1]) + wg_ij[-1, :] * np.sin(theta[-1])
)
a_ij[N, N] = np.sum(s_ij)

I = np.zeros(N + 1)
I[:N] = -(-np.cos(alpha) * np.sin(theta) + np.sin(alpha) * np.cos(theta))
I[N] = -(
    np.cos(alpha) * (np.cos(theta[0]) + np.cos(theta[-1])) +
    np.sin(alpha) * (np.sin(theta[0]) + np.sin(theta[-1]))
)

qg = np.linalg.solve(a_ij, I)

# plt.figure()
# plt.plot(x_pc, qg[:-1], "red")
# plt.xlim(0,1)
# plt.xticks(np.arange(0, 1.1, 0.1))
# plt.grid(True, linestyle='--')
# plt.legend()

b_ij = uq_ij * np.cos(theta[:, None]) + wq_ij * np.sin(theta[:, None])
d_ij = ug_ij * np.cos(theta[:, None]) + wg_ij * np.sin(theta[:, None])
b_ij = np.hstack((b_ij, np.sum(d_ij, axis=1, keepdims=True)))
t_i = np.cos(alpha) * np.cos(theta) + np.sin(alpha) * np.sin(theta)

vel = b_ij @ qg + t_i
cp = 1 - vel**2
cl = sum(-cp*cuerda*np.cos(theta - alpha))
cmc4 = sum((cp*(x_pc - 0.25)*np.cos(theta) + cp*y_pc*np.sin(theta)) * cuerda)
# cdp = sum(cp*cuerda*np.sin(theta - alpha))

# plt.figure()
# plt.plot(x_pc, qg[:-1])
# plt.xlabel("x_punto_control")
# plt.ylabel("qi")
# plt.xlim(0,1)
# plt.xticks(np.arange(0, 1.1, 0.1))
# plt.grid(True, linestyle='--')
# plt.legend()

plt.figure()
plt.plot(x[:N // 2], cp[:N // 2], "orange", label="Intradós")
plt.plot(x[N // 2:N], cp[N // 2:N], "blue", label="Extradós")
plt.gca().invert_yaxis()
plt.xlim(0,1)
# plt.ylim(1,-1.25)
plt.xlabel("x/c")
plt.ylabel("Cp")
plt.xticks(np.arange(0, 1.1, 0.1))
plt.grid(True, linestyle='--')
plt.legend()

print("Cl =",cl)
print("Cmc/4 =",cmc4)
print(qg[N])



