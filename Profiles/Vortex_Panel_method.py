# -*- coding: utf-8 -*-
"""
Created on Sun Mar 23 18:59:04 2025

@author: mdeto
"""

import numpy as np
import sympy as sp
import matplotlib.pyplot as plt

f, xf, t, N = 0.02, 0.3, 0.19, 200
x = sp.symbols('x')
zc1 = f / xf**2 * (2 * xf * x - x**2)
zc2 = f / (1 - xf)**2 * (1 - 2 * xf + 2 * xf * x - x**2)
zt = 5 * t * (0.2969 * sp.sqrt(x) + -0.126 * x -0.3516 * x**2 + 0.2843 * x**3 -0.1015 * x**4)
zc_ext = sp.Piecewise((zc1 + zt, x < xf), (zc2 + zt, x >= xf))
zc_int = sp.Piecewise((zc1 - zt, x < xf), (zc2 - zt, x >= xf))

beta = np.zeros(N + 1)
xn = 0.5 * (1 + np.cos(beta)) 
z = np.array([float(zc_int.subs(x, xn[0]))])  
x_pc, z_pc, theta, cuerda = np.zeros(N), np.zeros(N), np.zeros(N), np.zeros(N)

for i in range(N):
    beta[i + 1] = beta[i] + 2 * np.pi / N
    xn[i + 1] = 0.5 * (1 + np.cos(beta[i + 1]))  
    x_pc[i] = (xn[i + 1] + xn[i]) / 2 
    if beta[i + 1] >= np.pi:
        z = np.append(z, float(zc_ext.subs(x, xn[i + 1])))  
    else:
        z = np.append(z, float(zc_int.subs(x, xn[i + 1]))) 
    z_pc[i] = (z[i + 1] + z[i]) / 2
    theta[i] = np.arctan2(z[i + 1] - z[i], xn[i + 1] - xn[i])  
    if theta[i] < 0:
        theta[i] += 2 * np.pi
    cuerda[i] = (xn[i + 1] - xn[i]) * np.cos(theta[i]) + (z[i + 1] - z[i]) * np.sin(theta[i])

gamma = sp.symbols('gamma', real=True)
gamma = np.zeros(N)

x_ij, z_ij, r1, r2, theta1, theta2 = np.zeros((N, N)), np.zeros((N, N)), np.zeros((N, N)), np.zeros((N, N)), np.zeros((N, N)), np.zeros((N, N))
u_coef, w_coef, u_ij, w_ij = np.zeros((N, N)), np.zeros((N, N)), np.zeros((N, N)), np.zeros((N, N))

for i in range(N):
    for j in range(N):
        x_ij[i, j] = (x_pc[i] - xn[j]) * np.cos(theta[j]) + (z_pc[i] - z[j]) * np.sin(theta[j])
        z_ij[i, j] = -(x_pc[i] - xn[j]) * np.sin(theta[j]) + (z_pc[i] - z[j]) * np.cos(theta[j])
        r1[i, j] = np.sqrt(x_ij[i, j]**2 + z_ij[i, j]**2)
        r2[i, j] = np.sqrt((x_ij[i, j] - cuerda[j])**2 + z_ij[i, j]**2)
        theta1[i, j] = np.arctan2(z_ij[i, j], x_ij[i, j])
        theta2[i, j] = np.arctan2(z_ij[i, j], (x_ij[i, j] - cuerda[j])) 
        if i == j:
            theta1[i, j] = 0
            theta2[i, j] = np.pi

for i in range(N):
    for j in range(N):
        u_coef[i, j] = (theta2[i, j] - theta1[i, j]) / (2 * np.pi)
        if i == j:
            u_coef[i, j] = 0.5
        w_coef[i, j] = (1 / (2 * np.pi)) * np.log(r2[i, j] / r1[i, j])
        u_ij[i, j] = u_coef[i, j] * np.cos(theta[j]) - w_coef[i, j] * np.sin(theta[j])
        w_ij[i, j] = u_coef[i, j] * np.sin(theta[j]) + w_coef[i, j] * np.cos(theta[j])

alfa, v_inf = 0 * np.pi / 180, 1
I, t_i = np.zeros(N), np.zeros(N)

for i in range(N):
    I[i] = -v_inf * (-np.cos(alfa) * np.sin(theta[i]) + np.sin(alfa) * np.cos(theta[i]))
    t_i[i] = v_inf * (np.cos(alfa) * np.cos(theta[i]) + np.sin(alfa) * np.sin(theta[i]))

I[N - 1] = 0 

a_ij = np.zeros((N, N))
for i in range(N):
    for j in range(N):
        a_ij[i, j] = -u_ij[i, j] * np.sin(theta[i]) + w_ij[i, j] * np.cos(theta[i])

a_ij[N - 1, 0] = 1
a_ij[N - 1, N - 1] = 1
for j in range(1, N - 1):
    a_ij[N - 1, j] = 0

b_ij = np.zeros((N, N))
for i in range(N):
    for j in range(N):
        b_ij[i, j] = u_ij[i, j] * np.cos(theta[i]) + w_ij[i, j] * np.sin(theta[i])

gamma = np.linalg.inv(a_ij).dot(I)

vel = b_ij.dot(gamma) + t_i
vel_p, gamma_p, Cp, Cl = np.zeros(N), np.zeros(N), np.zeros(N), np.zeros(N)

for i in range(1, N):
    vel_p[i] = (vel[i] + vel[i - 1]) / 2
    gamma_p[i] = (gamma[i] + gamma[i - 1]) / 2
    Cp[i] = 1 - (vel_p[i] / v_inf)**2
    Cl[i] = 2 * gamma_p[i] * cuerda[i] / v_inf

vel_p[0] = (vel[0] + vel[N - 1]) / 2
gamma_p[0] = (gamma[0] + gamma[N - 1]) / 2

plt.figure(1)
plt.plot(x_pc, Cp)
plt.xlabel("x")
plt.ylabel("Cp")
plt.title("Pressure Coefficient vs x")
plt.gca().invert_yaxis()

# plt.figure(2)
# plt.plot(x_pc, vel_p)
# plt.xlabel("x")
# plt.ylabel("Velocity (m/s)")
# plt.title("Velocity Distribution")

Cl_t = np.sum(Cl)
# print("Total lift coefficient (Cl):", Cl_t)
print(theta1)
