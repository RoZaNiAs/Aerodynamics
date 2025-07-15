# -*- coding: utf-8 -*-
"""
Created on Wed Apr 30 16:24:36 2025

@author: mdeto
"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


f, xf = 0.02, 0.3
t = 0.19
cr, ct = 1, 0.5
b = 2
delta = 45


def zc(c, n=60):
    x = np.linspace(0, c, n)
    y = np.where(
        x <= xf * c,
        (f * c / (xf**2)) * (2 * xf * x - (x / c)**2 * c),
        (f * c / (1 - xf)**2) * ((1 - 2 * xf) + 2 * xf * x / c - (x / c)**2) * c
    )
    yt = 5 * t * (
        0.2969 * np.sqrt(x / c) - 0.1260 * (x / c)
        - 0.3516 * (x / c)**2 + 0.2843 * (x / c)**3
        - 0.1015 * (x / c)**4
    ) * c
    yext = y + yt
    yint = y - yt
    return x, yext, yint


num_span = 17    
num_chord = 20
y_half = np.linspace(0, b/2, num_span)
chords = np.linspace(cr, ct, num_span)
delta_rad = np.radians(delta)

X, Y, Z_upper, Z_lower = [], [], [], []

for i in range(num_span):
    c = chords[i]
    y_pos = y_half[i]
    x, zext, zint = zc(c, num_chord)
    x_offset = y_pos * np.tan(delta_rad)
    x_swept = x + x_offset
    y_line = np.full_like(x, y_pos)

    X.append(x_swept)
    Y.append(y_line)
    Z_upper.append(zext)
    Z_lower.append(zint)


X = np.array(X)
Y = np.array(Y)
Z_upper = np.array(Z_upper)
Z_lower = np.array(Z_lower)


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')


ax.plot_surface(X, Y, Z_upper, color='steelblue', alpha=0.7, rstride=1, cstride=1)
ax.plot_surface(X, Y, Z_lower, color='steelblue', alpha=0.7, rstride=1, cstride=1)
ax.plot_surface(X, -Y, Z_upper, color='steelblue', alpha=0.7, rstride=1, cstride=1)
ax.plot_surface(X, -Y, Z_lower, color='steelblue', alpha=0.7, rstride=1, cstride=1)


# ax.plot_wireframe(X, Y, Z_upper, color='black', linewidth=0.5)
# ax.plot_wireframe(X, Y, Z_lower, color='black', linewidth=0.5)
# ax.plot_wireframe(X, -Y, Z_upper, color='black', linewidth=0.5)
# ax.plot_wireframe(X, -Y, Z_lower, color='black', linewidth=0.5)

# Labels and appearance
ax.set_aspect("equal")
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")

plt.show()












