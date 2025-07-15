# -*- coding: utf-8 -*-
"""
Created on Thu Apr 10 01:11:41 2025

@author: mdeto
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def Fluent_reader(filename):
    x_ext, y_ext = [], []
    x_int, y_int = [], []

    with open(filename, 'r') as file:
        lines = file.readlines()[3:]  

    section = None

    for line in lines:
        line = line.strip()

        if line.startswith('((xy/key/label "extrados")'):
            section = 'extrados'
            continue
        elif line.startswith('((xy/key/label "intrados")'):
            section = 'intrados'
            continue
        elif line == ')':
            section = None
            continue

        if section in ['extrados', 'intrados']:
            if line and '(' not in line and ')' not in line:
                try:
                    parts = list(map(float, line.split()))
                    if section == 'extrados':
                        x_ext.append(parts[0])
                        y_ext.append(parts[1])
                    elif section == 'intrados':
                        x_int.append(parts[0])
                        y_int.append(parts[1])
                except ValueError:
                    pass  # Ignore lines that can't be parsed

    # Convert to NumPy arrays
    x_ext, y_ext = np.array(x_ext), np.array(y_ext)
    x_int, y_int = np.array(x_int), np.array(y_int)

    # Sort extrados
    sort_ext = np.argsort(x_ext)
    x_ext, y_ext = x_ext[sort_ext], y_ext[sort_ext]

    # Sort intrados
    sort_int = np.argsort(x_int)
    x_int, y_int = x_int[sort_int], y_int[sort_int]

    return x_ext, y_ext, x_int, y_int


def Xfoil_reader(x):
    data = np.loadtxt(x, skiprows=3)

    x = data[:, 0]
    y = data[:, 1]
    cp = data[:, 2]
    
    return x, y, cp

def Xfoil (x):
    cl = 0.1083*x + 0.2293
    cmc4 = 0.0007*x - 0.0424
    cd = 0.0004*x + 0.006
    return cl, cmc4, cd

xT = np.linspace(0, 10, 1000)

alpha = np.linspace(0, 10, 5)
cl = np.array([0.1967, 0.463, 0.7395, 0.9915, 1.2803])
cd = np.array([0.007, 0.007, 0.0078, 0.0090, 0.0123])
cm = np.array([-0.0367, -0.0338, -0.0333, -0.0290, -0.0321])
FluentCp_profile = "NACA_2319/ALFA_0/Fluent_data/cp0.txt"
XfoilCp_profile = "NACA_2319/ALFA_0/Xfoil_data/NACA_2319_0_ST.txt"
xFext, cpFext = Fluent_reader(FluentCp_profile)[0], Fluent_reader(FluentCp_profile)[1]
xFint, cpFint = Fluent_reader(FluentCp_profile)[2], Fluent_reader(FluentCp_profile)[3]
xX, yX, cpX = Xfoil_reader(XfoilCp_profile)

def cl_fit(x):
    coef = np.polyfit(alpha, cl, 1)
    return coef[0] * x + coef[1]

def cd_fit(x):
    coef = np.polyfit(alpha, cd, 1)
    return coef[0] * x + coef[1]

def cm_fit(x):
    coef = np.polyfit(alpha, cm, 1)
    return coef[0] * x + coef[1]

def clcd_fit(x):
    ratio = cl / cd
    coef = np.polyfit(alpha, ratio, 1)
    return coef[0] * x + coef[1]

xT = np.linspace(0, 10, 1000)

fig, ax1 = plt.subplots()
# plt.title("Cl, Cd, Cmc/4")

ax1.axhline(0, color='black', linewidth=1)
ax1.axvline(0, color='black', linewidth=1)
ax1.set_xlabel("α [º]")
ax1.set_ylabel("Cl")
ax1.set_ylim(-0.1, 1.6)
ax1.grid(True, which='both', axis='both')
plt.xlim(0, 10)

ax2 = ax1.twinx()
ax2.set_ylim(ax1.get_ylim()[0] * 100, ax1.get_ylim()[1] * 100)
ax2.set_ylabel("Cl/Cd")

line1 = ax1.plot(xT, cl_fit(xT), color='blue', label="Cl Fluent")[0]
line2 = ax1.plot(xT, cd_fit(xT), color='red', label="Cd Fluent")[0]
line3 = ax1.plot(xT, cm_fit(xT), color='orange', label="Cmc/4 Fluent")[0]
line4 = ax2.plot(xT, clcd_fit(xT), color='green', label="Cl/Cd Fluent")[0]

xfoil_cl = ax1.plot(xT, Xfoil(xT)[0], color="blue", linestyle="--", label="Cl XFOIL")[0]
xfoil_cd = ax1.plot(xT, Xfoil(xT)[2], color="red", linestyle="--", label="Cd XFOIL")[0]
xfoil_cm = ax1.plot(xT, Xfoil(xT)[1], color="orange", linestyle="--", label="Cmc/4 XFOIL")[0]
xfoil_clcd = ax2.plot(xT, Xfoil(xT)[0] / Xfoil(xT)[2], color="green", linestyle="--", label="Cl/Cd XFOIL")[0]

all_lines = [line1, xfoil_cl, line2, xfoil_cd, line3, xfoil_cm, line4, xfoil_clcd]
all_labels = [line.get_label() for line in all_lines]

ax1.legend(all_lines, all_labels, loc='upper left')


plt.tight_layout()
plt.show()


def Polar(x):
    return (np.polyfit(cl, cd, 2)[0])*x**2 + (np.polyfit(cl, cd, 2)[1])*x + (np.polyfit(cl, cd, 2)[2])

xPolar = np.linspace(0, 1.4, 1000)

# plt.figure()
# # plt.title("Polar")
# plt.scatter(cl, cd, color = "black", facecolor = "None", s = 75)
# plt.plot(xPolar, Polar(xPolar), "black")
# plt.xlim(0, 1.4)
# plt.grid(True)
# plt.xlabel("Cl")
# plt.ylabel("Cd")



# fig, ax1 = plt.subplots()

# ax1.plot(xFext, cpFext, color="red", label='Fluent')
# ax1.plot(xFint, cpFint, color="red")
# ax1.plot(xX, cpX, color='blue', linestyle='--', label='Xfoil')
# ax1.set_xlim(0, 1)
# ax1.set_ylim(1.25, -1)
# ax1.set_xlabel("x/c")
# ax1.set_ylabel("Cp")
# ax1.grid(True, which='both', axis='both')

# ax2 = ax1.twinx()

# y_offset = -1 
# y_scale = 2    
# airfoil_y_scaled = yX * y_scale + y_offset

# ax2.plot(xX, -airfoil_y_scaled, color='black')
# ax2.set_ylabel("y/c")
# ax2.set_ylim(ax1.get_ylim())
# yticks_cp = ax1.get_yticks()  
# yticks_airfoil = (-yticks_cp - y_offset) / y_scale  
# ax2.set_yticks(yticks_cp)
# ax2.set_yticklabels([f"{yt:.2f}" for yt in yticks_airfoil])
# ax1.legend()
# fig.tight_layout()
# plt.show()

# print("Fluent: cpmin =", min(cpFext), "en", xFext[np.argmin(cpFext)])
# print("XFOIL: cpmin =", min(cpX), "en", xX[np.argmin(cpX)])

# FluentWallStress0_profile = r"C:\Users\mdeto\Desktop\UNI\Tercero\Aerodinamica\Fluent\NACA_2319\ALFA_0\Fluent_data\xshearstress0.txt"
# FluentWallStress10_profile = r"C:\Users\mdeto\Desktop\UNI\Tercero\Aerodinamica\Fluent\NACA_2319\ALFA_10\Fluent_data\xshearstress10.txt"
# xsext0, sext0 = Fluent_reader(FluentWallStress0_profile)[0], Fluent_reader(FluentWallStress0_profile)[1]
# xsext10, sext10 = Fluent_reader(FluentWallStress10_profile)[0], Fluent_reader(FluentWallStress10_profile)[1]

# plt.figure()
# plt.plot(xsext0, sext0, label = "α = 0º")
# plt.plot(xsext10, sext10, label = "α = 10º")
# plt.xlabel("x/c")
# plt.ylabel("Esfuerzo en la pared [Pa]")
# plt.xlim(0 ,1)
# # plt.ylim(0,)
# plt.xticks(np.arange(0, 1.1, 0.1))
# plt.grid(True, linestyle='--')
# plt.legend()


