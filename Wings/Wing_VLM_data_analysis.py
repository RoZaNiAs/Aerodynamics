# -*- coding: utf-8 -*-
"""
Created on Sat Apr 26 12:26:34 2025

@author: mdeto
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

#NOTA DE USO
#ENTRE LAS LÍNEAS 244-278 HAY VARIAS FUNCIONES QUE EXPORTAN DATOS CALCULADOS A UN EXCEL. SE HA DISEÑADO DE ESTA FORMA PARA PODER VER DE FORMA MÁS VISUAL ALGUNAS MATRICES. 
#PARA PODER EMPLEARLAS CORRECTAMENTE SE REQUIERE CREAR UNA CARPETA EN EL DIRECTORIO LOCAL DONDE SE HA ALOJADO EL CÓDIGO LLAMADA TABLAS O ELIMINAR ESTA PARTE DEL DIRECTORIO EN EL CÓDIGO ENTRE LAS LÍNEAS DICHAS. 
#SI NO SE QUIERE GUARDAR NADA EN EXCEL Y SOLO COMPROBAR LAS GRÁFICAS, BASTA CON COMENTAR O ELIMINAR ESTAS LÍNEAS.

#DEFINICIÓN DE VARIABLES
A, S, lamda = 7, 1, 0.52
delta, epsilon = 2, -2.5
m3c4 = 0

deltar, epsilonr = delta * np.pi/180, epsilon * np.pi/180

b = np.sqrt(S*A)
cgm = S / b
cr = cgm / ( (1 + lamda) / 2)
ct = lamda * cr
cam = 2/3 * cr * ((1 + lamda + lamda**2) / (1 + lamda))
yca = (b/2)*((1 + 2*lamda) / (3*(1 + lamda)))
xca = 0.25*cr + np.tan(deltar)*yca

Nx = 1
Ny = 16

f = 0.02
xf = 0.3
t = 0.19


alpha = np.array([-12, -6, 0, 6, 12])
alpharad = alpha * (np.pi/180)

#FUNCIONES QUE PERMITEN PANELIZAR LAS ALAS

def PanelesDcha(N, alfa):
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
    
    bij = -(alfa + (panels3c4[:,1] * (epsilonr/ (b/2)) ) - m3c4)
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
    
    panels3c4, panelsc4, esquinas, esqtorb = np.zeros((N, 2)), np.zeros((N, 4)), np.zeros((N + 1, 6)), np.zeros((N, 4))
    esquinaizda, esquinadcha = np.zeros((N, 2)), np.zeros((N, 2))
    div = np.linspace(0, b/2, N + 1)
    esquinaizda, esquinadcha = ba(div), bs(div)
    esquinas[:,0], esquinas[:,1], esquinas[:,2] = esquinaizda, esquinadcha, div
    esquinas[:,3], esquinas[:,4], esquinas[:,5] = esquinaizda, esquinadcha, -div
    middle = 0.5*(div[1:] - div[:-1])[0]
    panels3c4[0][1] = middle
    esqtorb[:,0], esqtorb[:,1], esqtorb[:,2], esqtorb[:,3] = xc4(div)[:-1], div[:-1], xc4(div)[1:], div[1:]
    
    for i in range(N - 1):
        panels3c4[i + 1][1] = panels3c4[i][1] + 2*middle
        
    panels3c4[:,0] = x3c4(panels3c4[:,1])
    panelsc4[:,1], panelsc4[:,3] = div[:-1], div[1:]
    panelsc4[:,0], panelsc4[:,2] = xc4(panelsc4[:,1]), xc4(panelsc4[:,3])
    
    panelsc4_1, panelsc4_2 = np.zeros((N, 2)), np.zeros((N, 2))
    panelsc4_2[:,0], panelsc4_2[:,1] = panelsc4[:,0], -panelsc4[:,1]
    panelsc4_1[:,0], panelsc4_1[:,1] = panelsc4[:,2], -panelsc4[:,3]
    
    cij = np.zeros(N)
    for i in range(N):
        cij[i] = -0.5 * ((ba(div[i + 1]) - bs(div[i + 1])) + (ba(div[i]) - bs(div[i])))
    
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
    
    return V, panelsc4_1, panelsc4_2, panels3c4, bacoeff, bscoeff, xc4coeff, x3c4coeff, root, tip, xc4r, x3c4r, xc4t, x3c4t, esquinas, cij, middle, esqtorb

root, tip, xc4r, x3c4r, xc4t, x3c4t, esquinas, cij, middle, esqtorb = PanelesIzda(Ny)[8], PanelesIzda(Ny)[9], PanelesIzda(Ny)[10], PanelesIzda(Ny)[11], PanelesIzda(Ny)[12], PanelesIzda(Ny)[13], PanelesIzda(Ny)[14], PanelesIzda(Ny)[15], PanelesIzda(Ny)[16], PanelesIzda(Ny)[17]

#FUNCIÓN QUE PERMITEN OBTENER LOS COEFICIENTES AERODINÁMICOS

def Coefficients(N, alfa):
    cij = PanelesDcha(N, alfa)[2]
    middle = 2*PanelesDcha(N, alfa)[3]
    xc4 = PanelesDcha(N, alfa)[4][:,0]
    sij = PanelesDcha(N, alfa)[5]
    
    aij = PanelesDcha(N, alfa)[0] + PanelesIzda(N)[0]
    bij = PanelesDcha(N, alfa)[1]
    
    gamma = np.linalg.solve(aij, bij)
    
    cl = 2*gamma/cij
    cL = (2/S)*sum(cl*middle*cij)
    
    cmoyj = -cl* (xc4 /cij)
    cmoyw =  (2/ (S*cam)) * sum(cmoyj*sij*cij)
    cma = cmoyw + cL*(xca/cam)
    
    yc = PanelesDcha(N, alfa)[6][:,1]
    y1I = PanelesIzda(N)[1][:,1]
    y2D = PanelesDcha(N, alfa)[8][:,1]
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
    return cL, cdi, cmoyw, cma, cl, cdit, cmoyj

pc = PanelesIzda(Ny)[3]
cL, cdi, cmoyw, cma, cl, cdit, cmoyj = np.zeros(len(alpha)), np.zeros(len(alpha)), np.zeros(len(alpha)), np.zeros(len(alpha)), np.zeros((len(alpha), Ny)), np.zeros((len(alpha), Ny)), np.zeros((len(alpha), Ny))

for i in range(len(alpha)):
    cL[i] = Coefficients(Ny, alpharad[i])[0]
    cdi[i] = Coefficients(Ny, alpharad[i])[1]
    cmoyw[i] = Coefficients(Ny, alpharad[i])[2]
    cma[i] = Coefficients(Ny, alpharad[i])[3]
    cl[i] = Coefficients(Ny, alpharad[i])[4]
    cdit[i] = Coefficients(Ny, alpharad[i])[5]
    cmoyj[i] = Coefficients(Ny, alpharad[i])[6]

clalpha = np.polyfit(alpha, cL, 1)[0] * 180/np.pi
cl0 = np.polyfit(alpha, cL, 1)[1]
alpha0 = -cl0/clalpha * 180/np.pi

cliadbas = np.zeros((Ny, 2))
cliadbas[:,0] = (cl[0] - cl[1]) / (cL[0] - cL[1])
cliadbas[:,1] = cl[0] - cliadbas[:,0]*cL[0]

cdec = np.polyfit(alpha, cdi, 2)

bacoeff, bscoeff, xc4coeff, x3c4coeff = PanelesIzda(Ny)[4], PanelesIzda(Ny)[5], PanelesIzda(Ny)[6], PanelesIzda(Ny)[7]

cb = np.zeros((Ny,2))
cb[:,0], cb[:,1] = cij, np.full(Ny, middle)

#GUARDADO DE DATOS EN EXCEL. ESTAS FUNCIONES SE HAN REALIZADO CON AYUDA DE CHATGPT

df_cl = pd.DataFrame(cl, index=alpha, columns=[f'Panel {i+1}' for i in range(cl.shape[1])]) 
df_cl.index.name = 'Alpha (°)'
df_cl.to_excel("Tablas/Coeficientes_sustentacion_locales.xlsx")

df_cL = pd.DataFrame(cL, index=alpha, columns=["CL"])
df_cL.index.name = "Alpha"
df_cL.to_excel("Tablas/CL_vs_alpha.xlsx")

df_cliadbas = pd.DataFrame(cliadbas, columns=["Sustentación adicional", "Sustentación básica"])
df_cliadbas.index = [f"Panel {i+1}" for i in range(Ny)]
df_cliadbas.to_excel("Tablas/Coeficientes_sustentacion_adicional_basica.xlsx")

df_cdit = pd.DataFrame(cdit, index=alpha, columns=[f'Panel {i+1}' for i in range(cdit.shape[1])])
df_cdit.index.name = 'Alpha (°)'
df_cdit.to_excel("Tablas/Coeficientes_resistencia_locales.xlsx")

df_cdi = pd.DataFrame(cdi, index=alpha, columns=["CDi"])
df_cdi.index.name = "Alpha"
df_cdi.to_excel("Tablas/CDi_vs_alpha.xlsx")

columns = ["X izda", "X dcha", "Y", "X izda", "X dcha", "Y"]
df_esquina = pd.DataFrame(esquinas, columns=columns)
df_esquina.to_excel("Tablas/Esquinas_paneles.xlsx", index=False)

df_pc = pd.DataFrame(pc, columns=["X", "Y"])
df_pc.index = [f"Panel {i+1}" for i in range(Ny)]
df_pc.to_excel("Tablas/Puntos_de_control.xlsx")

df_esqtorb = pd.DataFrame(esqtorb, columns=["x1", "y1", "x2", "y2"])
df_esqtorb.index = [f"Torbellino {i+1}" for i in range(Ny)]
df_esqtorb.to_excel("Tablas/Esquinas_torbellinos.xlsx")

df_cuerdas = pd.DataFrame(cb, columns=["c_ij", "b_ij"])
df_cuerdas.index = [f"Panel {i+1}" for i in range(Ny)]
df_cuerdas.to_excel("Tablas/Cuerdas_panel.xlsx")

#EXPOSICIÓN DE LAS VARIABLES CALCULADAS

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
print("Datos de Panelización")
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
print("Resultados de Panelización")
print("-------------------------------------")
print("Perfil de raíz...")
print("Borde de ataque:", root[0])
print("Borde de salida:", root[1])
print("Punto c/4:", xc4r)
print("Punto 3c/4:", x3c4r)
print("Perfil de punta...")
print("Borde de ataque:", tip[0])
print("Borde de salida:", tip[1])
print("Punto c/4:", xc4t)
print("Punto 3c/4:", x3c4t)
print("Ángulo del borde de ataque (º) =", bacoeff[0] * 180/np.pi)
print("Ángulo en c/4 (º) =", xc4coeff[0] * 180/np.pi)
print("Ángulo en 3c/4 (º) =", x3c4coeff[0] * 180/np.pi)
print("Ángulo del borde de salida (º) =", bscoeff[0] * 180/np.pi)
print("-------------------------------------")
print("Resultados Aerodinámicos")
print("-------------------------------------")
print("Alpha(º) =", alpha)
print("CL =", cL)
print("Clalpha (rad-1) =", clalpha)
print("Clalpha=0 =", cl0)
print("AlphaCL=0 (º) =", alpha0)
print("CDi =", cdi)
print("Cd0 =",cdec[2])
print("Cd1 =",cdec[1])
print("Cd2 =",cdec[0])
print("CMW OY =", cmoyw)
print("CMW CA =", cma)
print("Pendiente CMW CA =", np.polyfit(alpha, cma, 1)[0])

NP = np.arange(1, Ny + 1, 1)
colors = ["purple", "red", "green", "orange", "blue"]

#CREACIÓN DE GRÁFICAS CON RESULTADOS

plt.figure()
for i in range(len(alpha)):
    plt.plot(pc[:, 1], cl[i], label=f'{int(alpha[i])}º', color=colors[i], marker="o")
    plt.plot(-pc[:, 1], cl[i], color=colors[i], marker="o")
plt.grid(True)
plt.xlabel("y [m]")
plt.ylabel("Cl(y)")
plt.xlim(-1.4, 1.4)
plt.axhline(0, color='black', linewidth=1) 
plt.axvline(0, color='black', linewidth=1)
plt.legend()
plt.show()

plt.figure()
plt.title("Distribución de coeficiente de sustentación total")
plt.plot(alpha, cL, marker = "o")
plt.grid("True")
plt.axhline(0, color='black', linewidth=1) 
plt.axvline(0, color='black', linewidth=1)
plt.xlabel("α(º)")
plt.ylabel("CL")

plt.figure()
plt.title("Distribución de coeficiente de sustentación adicional y básica")
plt.plot(pc[:,1], cliadbas[:,0], "blue", label = "Adicional", marker = "o")
plt.plot(-pc[:,1], cliadbas[:,0], "blue", marker = "o")
plt.plot(pc[:,1], cliadbas[:,1], "orange", label = "Básica", marker = "o")
plt.plot(-pc[:,1], cliadbas[:,1], "orange", marker = "o")
plt.grid("True")
plt.axhline(0, color='black', linewidth=1) 
plt.axvline(0, color='black', linewidth=1)
plt.xlabel("y [m]")
plt.ylabel("Cl(y)")
plt.xticks(np.linspace(-1.4, 1.4, 8))
plt.legend()

plt.figure()
plt.title("Distribución de coeficiente de resistencia local")
for i in range(len(alpha)):
    plt.plot(pc[:,1], cdit[i], color = colors[i], label=f'{int(alpha[i])}º', marker = "o")
    plt.plot(-pc[:,1], cdit[i], color = colors[i], marker = "o")
plt.grid("True")
plt.xlabel("y [m]")
plt.ylabel("CDi(y)")
plt.xticks(np.linspace(-1.4, 1.4, 8))
plt.axhline(0, color='black', linewidth=1) 
plt.axvline(0, color='black', linewidth=1)
plt.legend()

plt.figure()
plt.title("Distribución de coeficiente de resistencia total")
plt.plot(alpha, cdi, marker = "o")
plt.grid("True")
plt.axhline(0, color='black', linewidth=1) 
plt.axvline(0, color='black', linewidth=1)
plt.xlabel("α(º)")
plt.ylim(0,)
plt.ylabel("CDi")

plt.figure()
plt.title("Distribución de coeficiente de momento local")
for i in range(len(alpha)):
    plt.plot(pc[:,1], cmoyj[i], color = colors[i], label=f'{int(alpha[i])}º', marker = "o")
    plt.plot(-pc[:,1], cmoyj[i], color = colors[i], marker = "o")
plt.grid("True")
plt.xlabel("y [m]")
plt.ylabel("CMoY(y)")
plt.axhline(0, color='black', linewidth=1) 
plt.axvline(0, color='black', linewidth=1)
plt.xticks(np.linspace(-1.4, 1.4, 8))
plt.legend()

plt.figure()
plt.title("Distribución de coeficiente de momento respecto al eje y")
plt.plot(alpha, cmoyw, marker = "o")
plt.grid("True")
plt.axhline(0, color='black', linewidth=1) 
plt.axvline(0, color='black', linewidth=1)
plt.xlabel("α(º)")
plt.ylabel("CMoY")

plt.figure()
plt.title("Distribución de coeficiente de momento respecto al eje centro aerodinámico")
plt.plot(alpha, cma, marker = "o")
plt.grid("True")
plt.axhline(0, color='black', linewidth=1) 
plt.axvline(0, color='black', linewidth=1)
plt.xlabel("α(º)")
plt.ylabel("CMA")

plt.figure()
plt.title("Distribución de coeficiente de momento local respecto al eje y")
plt.plot(alpha, cmoyw, marker = "o")
plt.grid("True")
plt.axhline(0, color='black', linewidth=1) 
plt.axvline(0, color='black', linewidth=1)
plt.xlabel("α(º)")
plt.ylabel("CMoY")

plt.figure()
plt.title("Distribución de coeficiente de momento respecto al eje centro aerodinámico y al y")
plt.plot(alpha, cma, marker = "o", label = "CA")
plt.plot(alpha, cmoyw, marker = "o", label = "OY")
plt.grid("True")
plt.axhline(0, color='black', linewidth=1) 
plt.axvline(0, color='black', linewidth=1)
plt.xlabel("α(º)")
plt.ylabel("CM")
plt.legend()

plt.figure()
plt.title("Comparación entre XFLR5 y VLM (CL)")
plt.plot(alpha, cL, "black", label = "VLM")
plt.scatter(alpha, [-1.013, -0.546, -0.066, 0.415, 0.886], label="XFLR5 n=16", s=100, facecolors = "None", edgecolors='black')
plt.grid("True")
plt.axhline(0, color='black', linewidth=1) 
plt.axvline(0, color='black', linewidth=1)
plt.xlabel("α(º)")
plt.ylabel("CL")
plt.legend()

plt.figure()
plt.title("Comparación entre XFLR5 y VLM (CD)")
plt.plot(alpha, cdi, "black", label = "VLM")
plt.scatter(alpha, [0.045, 0.013, 0, 0.008, 0.035], label="XFLR5 n=16",s=100, facecolors = "None", edgecolors='black')
plt.grid("True")
plt.axhline(0, color='black', linewidth=1) 
plt.axvline(0, color='black', linewidth=1)
plt.xlabel("α(º)")
plt.ylabel("CL")
plt.legend()

plt.figure()
plt.title("Comparación entre XFLR5 y VLM (CMoY)")
plt.plot(alpha, cmoyw, "black",label = "VLM")
plt.scatter(alpha, [0.371, 0.202, 0.025, -0.152, -0.322], label="XFLR5 n=16",s=100, facecolors = "None", edgecolors='black')
plt.grid("True")
plt.axhline(0, color='black', linewidth=1) 
plt.axvline(0, color='black', linewidth=1)
plt.xlabel("α(º)")
plt.ylabel("CMoY")
plt.legend()
