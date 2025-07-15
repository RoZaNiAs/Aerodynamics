# -*- coding: utf-8 -*-
"""
Created on Sat Apr 26 12:27:25 2025

@author: mdeto
"""
import numpy as np
import matplotlib.pyplot as plt

A, S, lamda = 4, 1, 0.5
delta, epsilon = 45, -4
deltar, epsilonr = delta * np.pi/180, epsilon * np.pi/180

b = np.sqrt(S*A)
cgm = S / b
cr = cgm / ( (1 + lamda) / 2)
ct = lamda * cr
cam = 2/3 * cr * ((1 + lamda + lamda**2) / (1 + lamda))
yca = (b/2)*((1 + 2*lamda) / (3*(1 + lamda)))
xca = 0.25*cr + np.tan(deltar)*yca

Nx = 1
Ny = 2

f = 0.02
xf = 0.4
t = 0.19
m3c4 = 0

alpha = 5
alphar = alpha * (np.pi/180)

def PanelesDcha(N):
    root, tip = np.zeros((2, 2)), np.zeros((2, 2))
    root[:,0] = np.linspace(0, cr, 2)
    root[:,1] = 0
    tip[:,0] = np.linspace(0.25*cr + b/2 * np.tan(deltar) - 0.25*ct, 0.25*cr + b/2 * np.tan(deltar) + 0.75*ct, 2)
    tip[:,1] = b/2
    xc4r, x3c4r = 0.25*(root[:,0][1] - root[:,0][0]), 0.75*(root[:,0][1] - root[:,0][0])
    xc4t, x3c4t = tip[0][0] + 0.25*(tip[:,0][1] - tip[:,0][0]), tip[0][0] + 0.75*(tip[:,0][1] - tip[:,0][0])
    xc4coeff = np.polyfit([0, b/2], [xc4r, xc4t], 1)
    x3c4coeff = np.polyfit([0, b/2], [x3c4r, x3c4t], 1)
    bacoeff = np.polyfit([0, b/2], [root[0][0], tip[0][0]], 1)
    bscoeff = np.polyfit([0, b/2], [root[1][0], tip[1][0]], 1)
    
    def ba(y):
        return bacoeff[0]*y + bacoeff[1]
    def bs(y):
        return bscoeff[0]*y + bscoeff[1]
    def xc4(y):
        return xc4coeff[0]*y + xc4coeff[1]
    def x3c4(y):
        return x3c4coeff[0]*y + x3c4coeff[1]
    
    panels3c4, panelsc4, pointsc4 = np.zeros((N, 2)), np.zeros((N, 4)), np.zeros((N, 2))
    div = np.linspace(0, b/2, N + 1)
    middle = 0.5*(div[1:] - div[:-1])[0]
    panels3c4[0][1] = middle
    for i in range(N - 1):
        panels3c4[i + 1][1] = panels3c4[i][1] + 2*middle
    panels3c4[:,0] = x3c4(panels3c4[:,1])
    panelsc4[:,1], panelsc4[:,3] = div[:-1], div[1:]
    panelsc4[:,0], panelsc4[:,2] = xc4(panelsc4[:,1]), xc4(panelsc4[:,3])
    pointsc4[:,1] = panels3c4[:,1]
    pointsc4[:,0] = xc4(pointsc4[:,1])
    cij, sij = np.zeros(N), np.zeros(N)
    for i in range(N):
        cij[i] = -0.5 * ((ba(div[i + 1]) - bs(div[i + 1])) + (ba(div[i]) - bs(div[i])))
    sij = cij*2*middle
    
    panelsc4_1, panelsc4_2 = np.zeros((N, 2)), np.zeros((N, 2))
    panelsc4_1[:,0], panelsc4_1[:,1] = panelsc4[:,0], panelsc4[:,1]
    panelsc4_2[:,0], panelsc4_2[:,1] = panelsc4[:,2], panelsc4[:,3]
    
    a, bj, c, d, g, h = np.zeros((N, N)), np.zeros((N, N)), np.zeros((N, N)), np.zeros((N, N)), np.zeros((N, N)), np.zeros((N, N))
    
    for i in range(N):
        for j in range(N):
            a[i][j] = panels3c4[i][0] - panelsc4_1[j][0]
            bj[i][j] = panels3c4[i][1] - panelsc4_1[j][1]
            c[i][j] = panels3c4[i][0] - panelsc4_2[j][0]
            d[i][j] = panels3c4[i][1] - panelsc4_2[j][1]
    e = np.sqrt(a**2 + bj**2)
    f = np.sqrt(c**2 + d**2)
    g = a - c
    h = bj - d
    k = (g*a + h*bj)/e - (g*c + h*d)/f
    l = -(1 + a/e)/bj + (1 + c/f)/d
    V = ( 1 / (4*np.pi) )* ( (k / (a*d - c*bj) ) + l)
    
    bij = -(alphar + (panels3c4[:,1] * (epsilonr/ (b/2)) ) - m3c4)
    return V, bij, cij, middle, pointsc4, sij, panels3c4, panelsc4_1, panelsc4_2

def PanelesIzda(N):
    root, tip = np.zeros((2, 2)), np.zeros((2, 2))
    root[:,0] = np.linspace(0, cr, 2)
    root[:,1] = 0
    tip[:,0] = np.linspace(0.25*cr + b/2 * np.tan(deltar) - 0.25*ct, 0.25*cr + b/2 * np.tan(deltar) + 0.75*ct, 2)
    tip[:,1] = -b/2
    xc4r, x3c4r = 0.25*(root[:,0][1] - root[:,0][0]), 0.75*(root[:,0][1] - root[:,0][0])
    xc4t, x3c4t = tip[0][0] + 0.25*(tip[:,0][1] - tip[:,0][0]), tip[0][0] + 0.75*(tip[:,0][1] - tip[:,0][0])
    xc4coeff = np.polyfit([0, b/2], [xc4r, xc4t], 1)
    x3c4coeff = np.polyfit([0, b/2], [x3c4r, x3c4t], 1)
    bacoeff = np.polyfit([0, -b/2], [root[0][0], tip[0][0]], 1)
    bscoeff = np.polyfit([0, -b/2], [root[1][0], tip[1][0]], 1)
    
    def ba(y):
        return bacoeff[0]*y + bacoeff[1]
    def bs(y):
        return bscoeff[0]*y + bscoeff[1]
    def xc4(y):
        return xc4coeff[0]*y + xc4coeff[1]
    def x3c4(y):
        return x3c4coeff[0]*y + x3c4coeff[1]
    
    panels3c4, panelsc4 = np.zeros((N, 2)), np.zeros((N, 4))
    div = np.linspace(0, b/2, N + 1)
    middle = 0.5*(div[1:] - div[:-1])[0]
    panels3c4[0][1] = middle
    for i in range(N - 1):
        panels3c4[i + 1][1] = panels3c4[i][1] + 2*middle
    panels3c4[:,0] = x3c4(panels3c4[:,1])
    panelsc4[:,1], panelsc4[:,3] = div[:-1], div[1:]
    panelsc4[:,0], panelsc4[:,2] = xc4(panelsc4[:,1]), xc4(panelsc4[:,3])
    cij, sij = np.zeros(N), np.zeros(N)
    for i in range(N):
        cij[i] = -0.5 * ((ba(div[i + 1]) - bs(div[i + 1])) + (ba(div[i]) - bs(div[i])))
    sij = cij*2*middle
    
    panelsc4_1, panelsc4_2 = np.zeros((N, 2)), np.zeros((N, 2))
    panelsc4_2[:,0], panelsc4_2[:,1] = panelsc4[:,0], -panelsc4[:,1]
    panelsc4_1[:,0], panelsc4_1[:,1] = panelsc4[:,2], -panelsc4[:,3]
    
    a, bj, c, d, g, h = np.zeros((N, N)), np.zeros((N, N)), np.zeros((N, N)), np.zeros((N, N)), np.zeros((N, N)), np.zeros((N, N))
    
    for i in range(N):
        for j in range(N):
            a[i][j] = panels3c4[i][0] - panelsc4_1[j][0]
            bj[i][j] = panels3c4[i][1] - panelsc4_1[j][1]
            c[i][j] = panels3c4[i][0] - panelsc4_2[j][0]
            d[i][j] = panels3c4[i][1] - panelsc4_2[j][1]
    e = np.sqrt(a**2 + bj**2)
    f = np.sqrt(c**2 + d**2)
    g = a - c
    h = bj - d
    k = (g*a + h*bj)/e - (g*c + h*d)/f
    l = -(1 + a/e)/bj + (1 + c/f)/d
    V = ( 1 / (4*np.pi) )* ( (k / (a*d - c*bj) ) + l)
    
    return V, panelsc4_1, panelsc4_2

def Coefficients(N):
    cij = PanelesDcha(N)[2]
    middle = 2*PanelesDcha(N)[3]
    xc4 = PanelesDcha(N)[4][:,0]
    sij = PanelesDcha(N)[5]
    
    aij = PanelesDcha(N)[0] + PanelesIzda(N)[0]
    bij = PanelesDcha(N)[1]
    
    gamma = np.linalg.solve(aij, bij)
    
    cl = 2*gamma/cij
    cL = (2/S)*sum(cl*middle*cij)
    
    cmoyj = -cl* (xc4/cij)
    cmoyw =  (2/ (S*cam)) * sum(cmoyj*sij*cij)
    cma = cmoyw + cL*(xca/cam)
    
    yc = PanelesDcha(N)[6][:,1]
    y1I = PanelesIzda(N)[1][:,1]
    y2D = PanelesDcha(N)[8][:,1]
    wliI, wliD = np.zeros((N,N)), np.zeros((N,N))
    for i in range(N):
        for j in range(N):
            wliD[j][i] = (1 / (2*np.pi))*(1 / (yc[i] - y2D[j]))
            wliI[j][i] = - (1 / (2*np.pi))*(1 / (yc[i] - y1I[j]))
    wlibj = wliD + wliI
    
    gammaNT = np.zeros(N)
    
    for i in range(N - 1):
        gammaNT[i] = gamma[i] - gamma[i + 1]
    gammaNT[-1] = gamma[-1]
    
    alphainf = gammaNT@wlibj
    alpha0 = alphainf/2
    cdit = -cl*alpha0
    cdi = (2/S) * sum(cdit*middle*cij)
    return cL, cdi, cmoyw, cma

cL, cdi, cmoyw, cma = Coefficients(Ny)

print("-------------------------------------")
print("Geometría alar")
print("-------------------------------------")
print("Alargamiento =", A)
print("Superficie Alar =", S)
print("Estrechamiento =", lamda)
print("Flecha 1/4 (º) =", delta)
print("Torsion lineal (º) =", epsilon)
print("-------------------------------------")
print("Ala trapecial")
print("-------------------------------------")
print("b =", b)
print("b/2 =", b/2)
print("cr =", cr)
print("ct =", ct)
print("cgm =", cgm)
print("cam =", cam)
print("xca =", xca)
print("yca =", yca)
print("-------------------------------------")
print("Panelado")
print("-------------------------------------")
print("Nx =", Nx)
print("Nx =", Ny)
print("Total de paneles semienvergadura =", Nx * Ny)
print("-------------------------------------")
print("Datos del perfil")
print("-------------------------------------")
print("f(%)=", f*100)
print("xf(%)=", xf*10)
print("t(%)=", t*100)
print("Pendiente en 3c/4 =", m3c4)
print("-------------------------------------")
print("Resultados")
print("-------------------------------------")
print("Alpha(º) =", alpha)
print("CL =", cL)
print("CDi =", cdi)
print("CMW OY =", cmoyw)
print("CMW CA =", cma)
