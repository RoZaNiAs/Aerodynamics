# -*- coding: utf-8 -*-
"""
Created on Sat Mar  8 19:19:02 2025

@author: mdeto
"""
import matplotlib.pyplot as plt
import numpy as np
import sympy as sp
from scipy.integrate import quad
from scipy.optimize import fsolve

DNI = str(52046872)

for i in range(len(DNI)):
    V = int(DNI[3])
    W = int(DNI[4])
    X = int(DNI[5])
    Y = int(DNI[6])
    Z = int(DNI[7])
    
f = int((1 + (3*V)/8))/100
xf = (20 + 10*int((2*W/8)))/100
t = (12 + int((8*X)/9))/100
ca = (20 + 5*int((3*Y)/8))/100
deltaa = (5 + 5*int((2*Z)/8)) * (np.pi / 180)

#%% Funciones para graficar el perfil

def Perfil(x, timon):
    y = np.ones(len(x))
    if timon == True:
        for i in range(len(x)):
            if (x[i] <= xf):
                y[i] = (f / (xf**2)) * (2 * xf * x[i] - x[i]**2)
            if (xf < x[i] <= (1 - ca)):
                y[i] = (f / (1 - xf)**2) * ((1 - 2 * xf) + 2 * xf * x[i] - x[i]**2)
            if (x[i] > (1 - ca)):
                y[i] = ((f / (1 - xf)**2) * ((1 - 2 * xf) + 2 * xf * x[i] - x[i]**2)) - (x[i]- (1 - ca))*np.tan(deltaa)
    else:
        y = np.where(
            x <= xf,
            (f / (xf**2)) * (2 * xf * x - x**2),
            (f / (1 - xf)**2) * ((1 - 2 * xf) + 2 * xf * x - x**2)
        )

    yt = 5 * t * (0.2969 * np.sqrt(x) - 0.1260 * x - 0.3516 * x**2 + 0.2843 * x**3 - 0.1015 * x**4)
    yext = y + yt
    yint = y - yt
    
    return x, y, yext, yint

xP = np.linspace(0 ,1, 1000)
yLM = Perfil(xP, False)[1]
yext = Perfil(xP, False)[2]
yint = Perfil(xP, False)[3]

yLMT = Perfil(xP, True)[1]
yextT = Perfil(xP, True)[2]
yintT = Perfil(xP, True)[3]

plt.figure()
plt.plot(xP, yLM, color="red", linestyle = "-.", label="Línea Media")
plt.plot(xP, yext, color="black", label="Espesor")
plt.plot(xP, yint, color="black")
plt.plot(xP, yLMT, color="red", linestyle="-.")
plt.plot(xP, yextT, color="black", linestyle = "--")
plt.plot(xP, yintT, color="black", linestyle = "--")
plt.axhline(0, color='gray', linewidth=1)
plt.axvline(0, color='gray', linewidth=1)
# plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.xlim(0, 1)
plt.axis("scaled")
plt.xlabel("x/c")
plt.ylabel("y/c")
plt.xlim(0, 1)
plt.legend()
plt.show()

print("-------------------------------------------------------------------")
print("NACA",str(int(f*100)) + str(int(xf/10 * 100)) + str(int(t*100)))
print("-------------------------------------------------------------------")
print("f/c = ",f)
print("xf(%) = ",xf)
print("t(%) = ",t)
print("ca(%) = ",ca)
print("delta = ",deltaa / (np.pi / 180))

#%% Funciones de Theodorsen

def Theodorsen (x, alpha):
    alpha = alpha * (np.pi / 180)

    def epsilon(x):
        return A1 * np.sin(x - phi01) - A2 * np.sin(2 * x - phi02)

    def psi(x):
        return A1 * np.cos(x - phi01) - A2 * np.cos(2 * x - phi02) + psi0

    plt.figure()
    plt.plot(x, epsilon(x), 'green', label="Epsilon")
    plt.plot(x, psi(x), 'red', label="Psi")
    plt.axhline(0, color='gray', linewidth=1)
    plt.axvline(0, color='gray', linewidth=1)
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.xlim(0, 2 * np.pi)
    plt.legend()
    plt.show()

    a = 0.2452
    psi_vals = psi(x)
    eps_vals = epsilon(x)
    
    theta = x - eps_vals

    x_cart = 2 * a * np.cosh(psi_vals) * np.cos(theta)
    yp = 2 * a * np.sinh(psi_vals) * np.sin(theta)
    xp = -x_cart + 0.508

    plt.figure()
    plt.plot(xp, yp, 'black')
    plt.axhline(0, color='gray', linewidth=1)
    plt.axvline(0, color='gray', linewidth=1)
    plt.xticks(np.arange(0, 1.1, 0.1))
    plt.yticks(np.arange(-0.5, 0.55, 0.1))
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.xlabel("x/c")
    plt.ylabel("y/c")
    plt.ylim(-0.5, 0.5)
    plt.axis("equal")
    plt.xlim(0, 1)
    plt.show()

    def Depsilon(x):
        return A1 * np.cos(x - phi01) - 2 * A2 * np.cos(2 * x - phi02)

    def Dpsi(x):
        return -A1 * np.cos(x - phi01) - 2 * A2 * np.cos(2 * x - phi02)

    epsbs = epsilon(np.pi)
    Deps = Depsilon(x)
    Dps = Dpsi(x)

    epstheta = Deps * (1 + Deps / (1 - Deps))
    psitheta = Dps * (1 + Deps / (1 - Deps))

    C1 = np.sin(alpha + x) + np.sin(alpha + epsbs)
    C2 = ((1 + epstheta) * np.exp(psi0)) / (
        np.sqrt((np.sinh(psi_vals)**2 + np.sin(theta)**2) * (1 + psitheta**2))
    )

    cp = 1 - (C1 * C2)**2
    u = C1 * C2
    uextra = np.where(u > 0, u, 0)
    uintra = np.where(u <= 0, u, 0)

    return xp, epsbs, cp


A1 = abs(0.005*(t*100)+ 0.005)
A2 = abs(-0.00028*(t*100) - 0.00252)
psi0 = abs(0.00913*(t*100) - 0.00458)
phi01 = phi02 = 0
c = 1
phi = np.linspace (0, 360, 360) * (np.pi/180)
alpha = 0

plt.figure()
plt.plot(Theodorsen(phi, alpha)[0], Theodorsen(phi, alpha)[2], 'black')
plt.axhline(0, color='gray', linewidth=1)
plt.axvline(0, color='gray', linewidth=1)
plt.ylim(-0.8, 0.8)
plt.gca().invert_yaxis()
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.xlabel("x/c")
plt.ylabel("Cp")
plt.xlim(0, 1)
# plt.legend()
plt.show()

Clalpha = 2*np.pi*np.exp(psi0)
Cl = 2*np.pi*np.exp(psi0)*(alpha + Theodorsen(phi, 0)[1])
cpTh = Theodorsen(phi, 0)[2]

print("-------------------------------------------------------------------")
print("Theodorsen Values")
print("-------------------------------------------------------------------")
print("A1 =", A1)
print("A2 =", A2)
print("psi0 =", psi0)
print("Clalpha =", Clalpha)
print("Cl =", Cl)

#%% TPL
def TPL(x, alpha, n):
    alpha = alpha * (np.pi/180)
    def DLineaMedia1(x, n):
        theta = c / 2 * (1 + np.cos(x))
        return (f / xf**2) * (2 * xf - 2 * theta) * np.cos(n * x)

    def DLineaMedia2(x, n):
        theta = 0.5 * (1 + np.cos(x))
        return (f / (1 - xf)**2) * (2 * xf - 2 * theta) * np.cos(n * x)
    An = np.zeros(n + 1)
    An[0] = (-1 / np.pi) * quad(DLineaMedia2, 0, np.arccos(2 * xf - 1),args=(0))[0] + (-1 / np.pi) * quad(DLineaMedia1, np.arccos(2 * xf - 1), np.pi,args=(0))[0]
    for i in range(1, n + 1):
        An[i] = (-2 / np.pi) * quad(DLineaMedia2, 0, np.arccos(2 * xf - 1), args=(i))[0] + (-2 / np.pi) * quad(DLineaMedia1, np.arccos(2 * xf - 1), np.pi, args=(i,))[0]

    theta = np.arccos( 2*x- 1 )
    thetaf = np.arccos( 2*xf- 1 )
    clsum = np.zeros(len(x))
    for t in range(1, n):
        clsum += An[t] * np.sin(t * theta)
    cl = 4 *((alpha + An[0]) * np.tan(theta / 2) + clsum)
    cp = -cl/2
    
    return An, cp, thetaf

xTh = Theodorsen (phi, 0)[0]
n = 200

alphaid = -TPL (xTh, 0 , 0)[0] * (180/np.pi)
alpha = alphaid
A0 = TPL (xTh, 0 , 0)[0][0]
A1 = TPL (xTh, 0 , 1)[0][1]
A2 = TPL(xTh, 0, 2)[0][2]
An = TPL(xTh, 0, 200)[0]
cpTPL = TPL(xTh, alpha, n)[1]

plt.figure()
plt.plot(xTh, cpTPL, "black")
plt.axhline(0, color='gray', linewidth=1)  
plt.axvline(0, color='gray', linewidth=1)  
plt.xticks(np.arange(0, 2.1, 0.1)) 
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.xlabel("x/c")
plt.ylabel("Cp")
plt.gca().invert_yaxis()
plt.xlim(0,1)
plt.ylim(0,-0.25)

print("-------------------------------------------------------------------")
print("TPL")
print("-------------------------------------------------------------------")
print("Thetaf =", TPL(xTh, alpha, n)[2] * (180/np.pi))
print("A0 =", A0)
print("A1 =", A1)
print("A2 =", A2)
print("Cl0 =", 2* np.pi * (A0 + 0.5 * A1))
print("Ideal alpha(º) =", - A0 * (180/(np.pi)))
print("Ideal cl =", np.pi * A1 )
print("Cmca =", -(np.pi/4) * (A1 + A2))

#%% Theodorsen + TPL


def op1(x):
    return Clalpha*(np.pi/180)*x + (A0 + 0.5 * A1) * Clalpha

def op2(x):
    return 2* np.pi *(np.pi/180)*x + 2* np.pi * (A0 + 0.5 * A1)

def op3(x):
    return Clalpha*(np.pi/180)*x

plt.figure()
# plt.title(f"Cp Theodorsen + TPL Ideal alpha")
plt.plot(xTh, cpTPL + cpTh, 'blue', label="Extradós")
plt.plot(xTh, - cpTPL + cpTh, 'orange', label="Intradós")
plt.axhline(0, color='gray', linewidth=1)  
plt.axvline(0, color='gray', linewidth=1)  
plt.xticks(np.arange(0, 2.1, 0.1)) 
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.gca().invert_yaxis()
plt.xlabel("x/c")
plt.ylabel("Cp")
plt.xlim(0,1)
# plt.ylim(1,-1)
plt.legend()


plt.figure()
plt.plot(np.linspace(-2, 8, 1000), op1(np.linspace(-2, 8, 1000)), label = "Theodorsen + TPL")
plt.plot(np.linspace(-2, 8, 1000), op2(np.linspace(-2, 8, 1000)), label = "TPL")
plt.plot(np.linspace(-2, 8, 1000), op3(np.linspace(-2, 8, 1000)), label = "Theodorsen")
plt.axhline(0, color='black', linewidth=1)
plt.axvline(0, color='black', linewidth=1)
plt.xlabel("α[º]")
plt.ylabel("Cl")
plt.xlim(-int((A0 + 0.5 * A1)*(180/np.pi)+1),8)
plt.ylim(0,)
plt.legend()
plt.grid(True)
plt.show()

cpminTP = (min(cpTPL + cpTh))
xcpminTP = xTh[np.argmin(cpTPL + cpTh)]
xcenpTPL = (np.pi/4)*(A1 + A2) / (op1(np.linspace(-2, 6, 1000)))


print("-------------------------------------------------------------------")
print("Theodorsen + TPL")
print("-------------------------------------------------------------------")
print("Clalpha =", Clalpha)
print("Alpha0(º) =", -(A0 + 0.5 * A1) * 180/np.pi)
print("Cl0 =", (A0 + 0.5 * A1) * Clalpha)
print("Cmc/4 =", -(np.pi/4)*(A1 + A2))
print("xca = 0.25")
print("Cmca =", -(np.pi/4)*(A1 + A2))
print("Cpmin =",cpminTP,"en",xcpminTP)

#%%Paneles

def Paneles(alpha, timon, N):
    beta = np.linspace(0, 2 * np.pi, N + 1)
    x = 0.5 * (1 + np.cos(beta))
    alpha = alpha * (np.pi/180)

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
    
    return x, cp, cl, cmc4

N = 220
timon = False
alpha = alphaid

cpPan = Paneles(alpha, timon, N)[1]
cpminPan = min(Paneles(alpha, timon, N)[1])
xPan = Paneles(alphaid, timon, N)[0]
xcpminPan = Paneles(alpha, timon, N)[0][(np.argmin(Paneles(alpha, timon, N)[1]))]
clPan = Paneles(alpha, timon, N)[2]
cmc4Pan = Paneles(alpha, timon, N)[3]


plt.figure()
plt.plot(xPan[:N // 2], cpPan[:N // 2], "orange", label="Intradós")
plt.plot(xPan[N // 2:N], cpPan[N // 2:N], "blue", label="Extradós")
plt.gca().invert_yaxis()
plt.xlim(0,1)
plt.ylim(1,-1.25)
plt.xlabel("x/c")
plt.ylabel("Cp")
plt.xticks(np.arange(0, 1.1, 0.1))
plt.grid(True, linestyle='--')
plt.legend()
plt.show()

balpha = np.arange(-2, 7, 1) 
lotalpha = np.linspace(-2, 6, 100) 
clbalpha, cmbalpha = np.zeros(len(balpha)), np.zeros(len(balpha))
cllotalpha, cmlotalpha = np.zeros(len(lotalpha)), np.zeros(len(lotalpha))

for i in range(len(balpha)):
    clbalpha[i] = Paneles(balpha[i], timon, N)[2]
    cmbalpha[i] = Paneles(balpha[i], timon, N)[3]
for i in range(len(lotalpha)):
    cllotalpha[i] = Paneles(lotalpha[i], timon, N)[2]
    cmlotalpha[i] = Paneles(lotalpha[i], timon, N)[3]

Cl0Pan = np.polyfit(balpha, clbalpha, 1)[1]
ClalphaPan = np.polyfit(balpha, clbalpha, 1)[0]

def CambioPolo(x):
    cm1 = cmbalpha[0] + clbalpha[0]*(x - 0.25)
    cm2 = cmbalpha[1] + clbalpha[1]*(x - 0.25)
    return cm1 - cm2

xcaPan = float(fsolve(CambioPolo, 0.25))
cmcaPan = float(cmbalpha[0] + clbalpha[0]*(xcaPan - 0.25))
xcenpPan = - cmlotalpha/cllotalpha

plt.figure()
plt.plot(balpha, clbalpha)
plt.axhline(0, color='black', linewidth=1)
plt.axvline(0, color='black', linewidth=1)
plt.xlabel("α[º]")
plt.ylabel("Cl")
plt.xlim(-2,6)
plt.ylim(0,1)
plt.grid(True, linestyle='--')

plt.figure()
plt.plot(balpha, cmbalpha)
plt.axhline(0, color='black', linewidth=1)
plt.axvline(0, color='black', linewidth=1)
plt.xlabel("α[º]")
plt.ylabel("Cmc/4")
plt.xlim(-2,6)
plt.ylim(-0.06, 0)
plt.grid(True, linestyle='--')

plt.figure()
plt.plot(balpha, clbalpha, label = "Cl")
plt.plot(balpha, cmbalpha, label = "Cmc/4")
plt.axhline(0, color='black', linewidth=1)
plt.axvline(0, color='black', linewidth=1)
plt.xlim(-2,6)
plt.xlabel("α[º]")
plt.ylim(-0.2, 1)
plt.grid(True, linestyle='--')
plt.legend()

fig, ax1 = plt.subplots()
ax1.plot(lotalpha, cllotalpha, label="Cl")
ax1.plot(lotalpha, cmlotalpha*10, label = "cmc/4")
ax1.plot(lotalpha, np.full(len(lotalpha), -0.046)*10, label = "cmca")
ax1.set_xlim(-2, 6)
ax1.set_ylim(-0.8, 1)
ax1.set_xlabel("α [º]")
ax1.set_ylabel("Cl")
ax1.axhline(0, color='black', linewidth=1)
ax1.axvline(0, color='black', linewidth=1)
ax1.grid(True)
ax1.grid(True, which='both', axis='both') 
ax2 = ax1.twinx()
ax2.set_ylabel("Cmc/4")  
ax2.set_ylim(ax1.get_ylim()[0] / 10, ax1.get_ylim()[1] / 10)  
ticks = ax1.get_yticks()
ax2.set_yticks(ticks / 10)
ax1.grid(True, which='both', axis='both') 
ax1.legend(loc='upper left', bbox_to_anchor=(1.05, 1), borderaxespad=0.)
# plt.title("Lift, Moment Coefficients and Aerodinamic Center Moment vs Angle of Attack PANELES")
plt.show()


print("-------------------------------------------------------------------")
print("Paneles")
print("-------------------------------------------------------------------")
print("Cl =", clPan)
print("Cl0 =", Cl0Pan)
print("Clalpha[grados-1]=", ClalphaPan)
print("Clalpha[rad-1]=", ClalphaPan * (180/np.pi))
print("alpha0=", - Cl0Pan / ClalphaPan)
print("Cmc/4 =", np.polyfit(balpha, cmbalpha, 1)[1])
print("xca =", xcaPan)
print("cmca =", cmcaPan)
print("Cpmin =", cpminPan,"en",xcpminPan)

# %% XFOIL
perfil = "Airfoil_data/NACA_2319_0.56_ST.txt"
data = np.loadtxt(perfil, skiprows=3)

x = data[:, 0]
y = data[:, 1]
cp = data[:, 2]

def Xfoil (x):
    cl = 0.1257*x + 0.2466
    cmc4 = -0.0027*x - 0.048
    return cl, cmc4

xd = np.linspace (-10, 20, 1000)

plt.figure()
plt.plot(xd, Xfoil(xd)[0])
plt.axhline(0, color='black', linewidth=1)
plt.axvline(0, color='black', linewidth=1)
plt.xlabel("α [º]")
plt.ylabel("Cl")
plt.xlim(-10,20)
plt.ylim(-1.5,3)
plt.grid(True, linestyle='--')
# plt.legend()
plt.show()

plt.figure()
plt.plot(xd, Xfoil(xd)[1])
plt.axhline(0, color='black', linewidth=1)
plt.axvline(0, color='black', linewidth=1)
plt.xlabel("α [º]")
plt.ylabel("Cmc/4")
plt.xlim(-10,20)
plt.ylim(-0.3,0.3)
plt.grid(True, linestyle='--')
# plt.legend()
plt.show()

plt.figure()
plt.plot(xd, np.full(len(xd), -0.0427))
plt.axhline(0, color='black', linewidth=1)
plt.axvline(0, color='black', linewidth=1)
plt.xlabel("α [º]")
plt.ylabel("Cmca")
plt.xlim(-10,20)
plt.ylim(-0.3,0.3)
plt.grid(True, linestyle='--')
# plt.legend()
plt.show()

fig, ax1 = plt.subplots()
ax1.plot(xd, Xfoil(xd)[0], label="Cl")
ax1.plot(xd, Xfoil(xd)[1]*10, label = "cmc/4")
ax1.plot(xd, np.full(len(xd), -0.0427)*10, label = "cmca")
ax1.set_xlim(-2, 6)
ax1.set_ylim(-0.8, 1)
ax1.set_xlabel("α [º]")
ax1.set_ylabel("Cl")
ax1.axhline(0, color='black', linewidth=1)
ax1.axvline(0, color='black', linewidth=1)
ax1.grid(True)
ax1.grid(True, which='both', axis='both') 
ax2 = ax1.twinx()
ax2.set_ylabel("Cmc/4")  
ax2.set_ylim(ax1.get_ylim()[0] / 10, ax1.get_ylim()[1] / 10)  
ticks = ax1.get_yticks()
ax2.set_yticks(ticks / 10)
ax1.grid(True, which='both', axis='both') 
ax1.legend(loc='upper left', bbox_to_anchor=(1.05, 1), borderaxespad=0.)
# plt.title("Lift, Moment Coefficients and Aerodinamic Center Moment vs Angle of Attack XFOIL")
plt.show()

#%%Paneles vs XFOIL vs TPL

perfil = "Airfoil_data/NACA_2319_0.56_ST.txt"
data = np.loadtxt(perfil, skiprows=3)

xX = data[:, 0]
yX = data[:, 1]
cpXFOIL = data[:, 2]

fig, ax1 = plt.subplots()
alpha_range = np.linspace(-2, 6, 1000)
ax1.plot(alpha_range, op1(alpha_range), label="Cl Theodorsen + TPL")
ax1.plot(lotalpha, cllotalpha, label="Cl Paneles")
ax1.plot(alpha_range, Xfoil(alpha_range)[0], label="Cl Xfoil")
ax1.plot(alpha_range, xcenpTPL + 0.25, "--", label="xcp TPL + Theodorsen")
ax1.plot(lotalpha, xcenpPan + 0.25, "--", color = "black", label="xcp Paneles")
ax1.plot(alpha_range, - ( Xfoil(alpha_range)[1] / Xfoil(alpha_range)[0]) + 0.25, "--", label="xcp Xfoil")
cmc_theodorsen = np.full_like(alpha_range, -(np.pi/4)*(A1 + A2))
ax1.plot(alpha_range, cmc_theodorsen*10, "-.", label="Cmc/4 Theodorsen + TPL")
ax1.plot(lotalpha, cmlotalpha*10, "-.", label="Cmc/4 Paneles")
ax1.plot(alpha_range, Xfoil(alpha_range)[1]*10, "-.", label="Cmc/4 Xfoil")
ax1.set_xlim(-2, 6)
ax1.set_ylim(-0.8, 1)
ax1.set_xlabel("α [º]")
ax1.set_ylabel("Cl")
ax1.axhline(0, color='black', linewidth=1)
ax1.axvline(0, color='black', linewidth=1)
ax1.grid(True)
ax1.grid(True, which='both', axis='both') 
ax2 = ax1.twinx()
ax2.set_ylabel("Cmc/4")  
ax2.set_ylim(ax1.get_ylim()[0] / 10, ax1.get_ylim()[1] / 10)  
ticks = ax1.get_yticks()
ax2.set_yticks(ticks / 10)
ax1.grid(True, which='both', axis='both') 
ax1.legend(loc='upper left', bbox_to_anchor=(1.05, 1), borderaxespad=0.)
# plt.title("Lift, Moment Coefficients and Pressure Center vs Angle of Attack")
plt.show()


plt.figure()
plt.scatter(xX, cpXFOIL, color="black", label = "XFOIL")
plt.plot(Paneles(alphaid, False, 220)[0][:-1], Paneles(alphaid, False, 220)[1],"red",label = "220 Paneles")
plt.plot(xTh, cpTPL + cpTh, 'blue', label="TPL")
plt.plot(xTh, - cpTPL + cpTh, 'blue')
plt.gca().invert_yaxis()
plt.xlim(0,1)
plt.ylim(1,)
plt.xlabel("x/c")
plt.ylabel("Cp")
plt.xticks(np.arange(0, 1.1, 0.1))
plt.grid(True, linestyle='--')
plt.legend()
plt.show()

print("XFOIL: cpmin =",min(cpXFOIL), "en", xX[np.argmin(cpXFOIL)])
#%% TPL Timon

def GlauertTimon(x, alpha, n):
    alpha = alpha * (np.pi/180)
    Ant = np.zeros(n + 1)
    thetaa = np.arccos(2 * (1-ca) - 1)
    def Dtimon(theta, n):
        return  - deltaa*np.cos(n*theta)
    Ant[0] = (-1 / np.pi) * (quad(Dtimon, 0, thetaa,args=(0))[0])
    for i in range(1, n + 1):
        Ant[i] = (-2 / np.pi) * (quad(Dtimon, 0, thetaa,args=i)[0])
    theta = np.arccos( 2*x- 1 )
    clsum = np.zeros(len(x))
    A = Ant + An
    A0 = A[0]
    for t in range(1, n):
        clsum += A[t] * np.sin(t * theta)
    cl = 4*((alpha + A0) * np.tan(theta / 2) + clsum)
    cp = -cl/2
    
    return cp, Ant, thetaa

GlauertTimon(xTh, 0, 200)
cpTimon = GlauertTimon(xTh, 0, 200)[0]
Ant = GlauertTimon(xTh, 0, 200)[1]
A = An + Ant

plt.figure()
# plt.title("Cp (alpha 0)")
plt.plot(xTh, cpTimon, "black", label = "200")
plt.axhline(0, color='gray', linewidth=1)  
plt.axvline(0, color='gray', linewidth=1)  
plt.xticks(np.arange(0, 2.1, 0.1)) 
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.xlabel("x/c")
plt.ylabel("Cp")
plt.gca().invert_yaxis()
plt.axis("scaled")
plt.xlim(0,1)
plt.ylim(0,-0.8)
plt.legend()
plt.show()

def ch(alpha, n):
    theta_a = np.arccos(2 * (1 - ca) - 1)
    A0 = A[0]
    def integrand(theta):
        clsum = sum(A[t] * np.sin(t * theta) for t in range(1, n))
        cl_theta = 4 * (alpha + A0) * np.tan(theta / 2) + 4 * clsum
        return (cl_theta / 4) * (np.cos(theta) - np.cos(theta_a)) * np.sin(theta)
    ch_value = quad(integrand, 0, theta_a)[0]
    return ch_value

ch = ch(0, 200)

plt.figure()
# plt.title(f"Cp Theodorsen + TPL Alpha 0 Timon")
plt.plot(xTh, cpTh + cpTimon, 'blue', label="Extradós")
plt.plot(xTh, cpTh - cpTimon, 'orange', label="Intradós")
plt.axhline(0, color='gray', linewidth=1)  
plt.axvline(0, color='gray', linewidth=1)  
plt.xticks(np.arange(0, 2.1, 0.1)) 
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.gca().invert_yaxis()
plt.xlabel("x/c")
plt.ylabel("Cp")
plt.xlim(0,1)
plt.ylim(1,-1.2)
plt.legend()


print("-------------------------------------------------------------------")
print("TPL Timón")
print("-------------------------------------------------------------------")
print("Thetaa =", GlauertTimon(xTh, 0, 200)[2] * (180/np.pi))
print("A0 Timón =", Ant[0])
print("A1 Timón =", Ant[1])
print("A2 Timón =", Ant[2])
print("A0 Total =", A[0])
print("A1 Total =", A[1])
print("A2 Total =", A[2])
print("Cl timón(alpha = 0º) =", 2* np.pi * (Ant[0] + 0.5 * Ant[1]))
print("Cmca timón(alpha = 0º) =", - np.pi / 4 * (Ant[1] + 0.5 * Ant[2]))
print("Cl completo(alpha = 0º) =", 2* np.pi * (A[0] + 0.5 * A[1]))
print("Cmca completo(alpha = 0º) =", - np.pi / 4 * (A[1] + (A[2])))
print("DeltaCl =", 2* np.pi * ( (A[0] + 0.5 * A[1]) - (An[0] + 0.5 * An[1])  ) )
print("DeltaCmca =", - (np.pi / 4) * ( A[1] + A[2] - (An[1] + An[2]) ) )
print("Ch =", ch)

#%% XFOIL vs TPL Timon

def ClTplTimon(x):
    return 2* np.pi *(np.pi/180)*x + 2* np.pi * (Ant[0] + 0.5 * Ant[1])

def ClTotalTimon(x):
    return Clalpha*(np.pi/180)*x + (A[0] + 0.5 * A[1]) * Clalpha

def ClTotalTimonXfoil(x):
    return 0.1247*x + 0.6727

def CmTotalTimonXfoil(x):
    return -0.0025*x - 0.1136

xTimon = np.linspace (-2, 6 , 1000)

plt.figure()
plt.title("Cl solo timon")
plt.plot(xTimon, ClTplTimon(xTimon))
plt.axhline(0, color='gray', linewidth=1)  
plt.axvline(0, color='gray', linewidth=1) 
plt.xlim(-2, 6)
plt.ylim(0 ,1)
plt.grid(True, which='both', linestyle='--', linewidth=0.5)

plt.figure()
plt.title("Cl solo timon")
plt.plot(xTimon, ClTplTimon(xTimon))
plt.axhline(0, color='gray', linewidth=1)  
plt.axvline(0, color='gray', linewidth=1) 
plt.xlim(-2, 6)
plt.ylim(0 ,1)
plt.grid(True, which='both', linestyle='--', linewidth=0.5)

plt.figure()
plt.title("Cmc/4 solo timon")
plt.plot(xTimon, np.full(len(xTimon), - np.pi / 4 * (Ant[1] + 0.5 * Ant[2])))
plt.axhline(0, color='gray', linewidth=1)  
plt.axvline(0, color='gray', linewidth=1) 
plt.xlim(-2, 6)
# plt.ylim(0 ,1)
plt.grid(True, which='both', linestyle='--', linewidth=0.5)

fig, ax1 = plt.subplots()
ax1.plot(xTimon, ClTplTimon(xTimon), label="Cl")
ax1.plot(xTimon, np.full(len(xTimon), - np.pi / 4 * (Ant[1] + 0.5 * Ant[2]))*10, label = "cmc/4")
ax1.set_xlim(-2, 6)
ax1.set_ylim(-0.8, 1)
ax1.set_xlabel("α [º]")
ax1.set_ylabel("Cl")
ax1.axhline(0, color='black', linewidth=1)
ax1.axvline(0, color='black', linewidth=1)
ax1.grid(True)
ax1.grid(True, which='both', axis='both') 
ax2 = ax1.twinx()
ax2.set_ylabel("Cmc/4")  
ax2.set_ylim(ax1.get_ylim()[0] / 10, ax1.get_ylim()[1] / 10)  
ticks = ax1.get_yticks()
ax2.set_yticks(ticks / 10)
ax1.grid(True, which='both', axis='both') 
ax1.legend(loc='upper left', bbox_to_anchor=(1.05, 1), borderaxespad=0.)
# plt.title("Lift, Moment Coefficients  vs Angle of Attack SOLO TIMON")
plt.show()

plt.figure()
plt.title("Cl total timon")
plt.plot(xTimon, ClTotalTimon(xTimon))
plt.axhline(0, color='gray', linewidth=1)  
plt.axvline(0, color='gray', linewidth=1) 
plt.xlim(-2, 6)
plt.ylim(0 ,1)
plt.grid(True, which='both', linestyle='--', linewidth=0.5)

plt.figure()
plt.title("Cmc/4 total timon")
plt.plot(xTimon, np.full(len(xTimon), - np.pi / 4 * (A[1] + 0.5 * A[2])))
plt.axhline(0, color='gray', linewidth=1)  
plt.axvline(0, color='gray', linewidth=1) 
plt.xlim(-2, 6)
# plt.ylim(0 ,1)
plt.grid(True, which='both', linestyle='--', linewidth=0.5)

fig, ax1 = plt.subplots()
ax1.plot(xTimon, ClTotalTimon(xTimon), label="Cl")
ax1.plot(xTimon, np.full(len(xTimon), - np.pi / 4 * (A[1] + 0.5 * A[2]))*10, label = "cmc/4")
ax1.set_xlim(-2, 6)
ax1.set_ylim(-1.5, 1.5)
ax1.set_xlabel("α [º]")
ax1.set_ylabel("Cl")
ax1.axhline(0, color='black', linewidth=1)
ax1.axvline(0, color='black', linewidth=1)
ax1.grid(True)
ax1.grid(True, which='both', axis='both') 
ax2 = ax1.twinx()
ax2.set_ylabel("Cmc/4")  
ax2.set_ylim(ax1.get_ylim()[0] / 10, ax1.get_ylim()[1] / 10)  
ticks = ax1.get_yticks()
ax2.set_yticks(ticks / 10)
ax1.grid(True, which='both', axis='both') 
ax1.legend(loc='upper left', bbox_to_anchor=(1.05, 1), borderaxespad=0.)
# plt.title("Lift, Moment Coefficients  vs Angle of Attack TOTAL TIMON")
plt.show()

plt.figure()
# plt.title("Cl total timon Xfoil")
plt.plot(xTimon, ClTotalTimonXfoil(xTimon))
plt.axhline(0, color='gray', linewidth=1)  
plt.axvline(0, color='gray', linewidth=1) 
plt.xlabel("α [º]")
plt.ylabel("Cl")
plt.xlim(-2, 6)
plt.ylim(0 ,1.6)
plt.grid(True, which='both', linestyle='--', linewidth=0.5)

plt.figure()
# plt.title("Cmc/4 total timon Xfoil")
plt.plot(xTimon, CmTotalTimonXfoil(xTimon))
plt.axhline(0, color='gray', linewidth=1)  
plt.axvline(0, color='gray', linewidth=1) 
plt.xlabel("α [º]")
plt.ylabel("Cmc/4")
plt.xlim(-2, 6)
plt.ylim(-0.3 ,0.3)
plt.grid(True, which='both', linestyle='--', linewidth=0.5)

fig, ax1 = plt.subplots()
ax1.plot(xTimon, ClTotalTimonXfoil(xTimon), label="Cl")
ax1.plot(xTimon, CmTotalTimonXfoil(xTimon)*10, label = "cmc/4")
ax1.set_xlim(-2, 6)
ax1.set_ylim(-1.5, 1.5)
ax1.set_xlabel("α [º]")
ax1.set_ylabel("Cl")
ax1.axhline(0, color='black', linewidth=1)
ax1.axvline(0, color='black', linewidth=1)
ax1.grid(True)
ax1.grid(True, which='both', axis='both') 
ax2 = ax1.twinx()
ax2.set_ylabel("Cmc/4")  
ax2.set_ylim(ax1.get_ylim()[0] / 10, ax1.get_ylim()[1] / 10)  
ticks = ax1.get_yticks()
ax2.set_yticks(ticks / 10)
ax1.grid(True, which='both', axis='both') 
ax1.legend(loc='upper left', bbox_to_anchor=(1.05, 1), borderaxespad=0.)
# plt.title("Lift, Moment Coefficients  vs Angle of Attack TOTAL TIMON XFOIL")
plt.show()

fig, ax1 = plt.subplots()
ax1.plot(xTimon, ClTotalTimon(xTimon), label="Cl TPL + Theodorsen")
ax1.plot(xTimon, np.full(len(xTimon), - np.pi / 4 * (A[1] + 0.5 * A[2]))*10, label = "cmc/4 TPL + Theodorsen")
ax1.plot(xTimon, ClTotalTimonXfoil(xTimon), label="Cl XFOIL")
ax1.plot(xTimon, CmTotalTimonXfoil(xTimon)*10, label = "cmc/4 XFOIL")
ax1.set_xlim(-2, 6)
ax1.set_ylim(-1.5, 1.5)
ax1.set_xlabel("α [º]")
ax1.set_ylabel("Cl")
ax1.axhline(0, color='black', linewidth=1)
ax1.axvline(0, color='black', linewidth=1)
ax1.grid(True)
ax1.grid(True, which='both', axis='both') 
ax2 = ax1.twinx()
ax2.set_ylabel("Cmc/4")  
ax2.set_ylim(ax1.get_ylim()[0] / 10, ax1.get_ylim()[1] / 10)  
ticks = ax1.get_yticks()
ax2.set_yticks(ticks / 10)
ax1.grid(True, which='both', axis='both') 
ax1.legend(loc='upper left', bbox_to_anchor=(1.05, 1), borderaxespad=0.)
# plt.title("Lift, Moment Coefficients  vs Angle of Attack TOTAL TIMON XFOIL")
plt.show()

perfil = "Airfoil_data/NACA_2319_0_T.txt"
data = np.loadtxt(perfil, skiprows=3)

xXt = data[:, 0]
yXt = data[:, 1]
cpXFOILt = data[:, 2]
timon = True
N = 220
cpPanT = Paneles(0, timon, N)[1]
xPanT = Paneles(0, timon, N)[0]

plt.figure()
# plt.title("TPL vs XFOIL vs PANELES TIMON")
plt.plot(xTh, cpTh + cpTimon, 'blue',  label= "TPL + Theodorsen")
plt.plot(xTh, cpTh - cpTimon, 'blue')
plt.scatter(xXt, cpXFOILt, color = "black", label = "XFOIL")
plt.plot(xPanT[:-1], cpPanT, "red", label = "220 Paneles" )
plt.axhline(0, color='gray', linewidth=1)  
plt.axvline(0, color='gray', linewidth=1)  
plt.xticks(np.arange(0, 2.1, 0.1)) 
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.gca().invert_yaxis()
plt.xlabel("x/c")
plt.ylabel("Cp")
plt.xlim(0,1)
plt.ylim(1,-1.25)
plt.legend()
