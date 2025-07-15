# -*- coding: utf-8 -*-
"""
Created on Mon Mar  3 16:45:28 2025

@author: mdeto
"""
import matplotlib.pyplot as plt
import numpy as np

def epsilon (x):
    return A1*np.sin(x - phi01) - A2*np.sin(2*x - phi02)

def Depsilon (x):
    return A1*np.cos(x - phi01) - 2*A2*np.cos(2*x - phi02)

def psi (x):
    return A1*np.cos(x - phi01) - A2*np.cos(2*x - phi02) + psi0

def Dpsi (x):
    return -A1*np.cos(x - phi01) - 2*A2*np.cos(2*x - phi02) 

    
def Theodorsen ():
    c=1
    theta = phi - epsilon(phi)
    x = (c/2) * np.cosh(psi(phi)) * np.cos(theta)
    y = (c/2) * np.sinh(psi(phi)) * np.sin(theta)
    return theta , epsilon, psi,  -x + 0.5, y


def cp(alpha):
    theta = phi - epsilon(phi)
    epsbs = epsilon(np.pi)
    epstheta = Depsilon(phi) * (1 + Depsilon(phi) / (1 - Depsilon(phi)))
    psitheta = Dpsi(phi) * (1 + Dpsi(phi) / (1 - Dpsi(phi)))
    C1 = np.sin(alpha + phi) + np.sin(alpha + epsbs)
    C2 = ( ( (1 + epstheta)*np.exp(psi0) ) / ( np.sqrt((np.sinh(psi(phi))**2 + np.sin(theta)**2) * (1 + psitheta**2) ) ) )
    cp = 1 - (C1*C2)**2
    return cp

def t():
    t = np.zeros(180)
    for i in range (179):
        t[i] = Theodorsen()[4][i] - Theodorsen()[4][180 + i]
    return max(t), Theodorsen()[3][np.argmax(t)]
    
    
    
A1 = 0.0665
A2 = abs(-0.004)
phi01 = 0
phi02 = 0
psi0 = 0.106
    
phi = np.linspace (0,360,360) * (np.pi/180)
alpha = 0 * (np.pi / 180)

Clalpha = 2*np.pi*np.exp(psi0)
Cl = 2*np.pi*np.exp(psi0)*(alpha + epsilon(np.pi))

print("Clalpha =", Clalpha)
print("Cl =", Cl)
print("Max t =", t()[0],"at x  =",t()[1])


plt.figure()
plt.plot(Theodorsen()[3], Theodorsen()[4], 'black')
plt.ylim(-0.2, 0.2)
plt.axhline(0, color='gray', linewidth=1)  
plt.axvline(0, color='gray', linewidth=1)  
plt.grid(True, which='both', linestyle='--', linewidth=0.5)


plt.figure()
plt.plot(Theodorsen()[3], cp(alpha), 'black')
plt.ylim(-1, 1)
plt.axhline(0, color='gray', linewidth=1)  
plt.axvline(0, color='gray', linewidth=1)  
plt.gca().invert_yaxis()
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
