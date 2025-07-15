# -*- coding: utf-8 -*-
"""
Created on Mon Apr  7 22:27:02 2025

@author: mdeto
"""
import numpy as np
import matplotlib.pyplot as plt

f, xf, t = 0.01, 0.2, 0.18

x = np.linspace(0, 1, 200)
z = np.where(
    x <= xf,
    (f / (xf**2)) * (2 * xf * x - x**2),
    (f / (1 - xf)**2) * ((1 - 2 * xf) + 2 * xf * x - x**2)
)
zt = 5 * t * (0.2969 * np.sqrt(x) - 0.1260 * x - 0.3516 * (x)**2 + 0.2843 * (x)**3 - 0.1015 * (x)**4)
zext = z + zt
zint = z - zt
zbs = np.array([zext[-1], zint[-1]])

puntos = np.linspace(1, len(x), len(x), dtype=int)

data1 = np.column_stack((np.ones(len(zext)), np.linspace(1, len(x), len(x), dtype=int), x, zext, np.zeros(len(zext))))
data2 = np.column_stack((np.full(len(zint), 2), np.linspace(1, len(x), len(x), dtype=int), x, zint, np.zeros(len(zint))))
data3 = np.column_stack((np.full(len(zbs), 3), np.linspace(1, len(zbs), len(zbs), dtype=int), np.ones(2), zbs, np.zeros(len(zbs))))

# Header without additional #
header = f"NACA_{f*100:.0f}{xf*10:.0f}{t*100:.0f}\nExtrados\nCurva        Punto      Coordenada x  Coordenada y  Coordenada z"

np.savetxt(
    f"NACA_{f*100:.0f}{xf*10:.0f}{t*100:.0f}.txt",
    data1,
    header=header,
    fmt="%10.6f   %10d   %10.6f   %10.6f   %10.6f"
)

with open(f"NACA_{f*100:.0f}{xf*10:.0f}{t*100:.0f}.txt", "a") as file:
    file.write("\n#Intrados\n") 
    np.savetxt(
        file,
        data2,
        fmt="%10.6f   %10d   %10.6f   %10.6f   %10.6f"
    )

with open(f"NACA_{f*100:.0f}{xf*10:.0f}{t*100:.0f}.txt", "a") as file:
    file.write("\n#Borde de salida\n")  
    np.savetxt(
        file,
        data3,
        fmt="%10.6f   %10d   %10.6f   %10.6f   %10.6f"
    )



plt.figure()
plt.plot(x,zext)
plt.plot(x,zint)
plt.plot(np.ones(2),zbs)
plt.axis("scaled")