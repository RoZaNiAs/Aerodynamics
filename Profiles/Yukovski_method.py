# -*- coding: utf-8 -*-
"""
Created on Sun Mar  2 15:04:50 2025

@author: mdeto
"""
import matplotlib.pyplot as plt
import numpy as np

def Yukowski(a, x, y):
    xn = x * (1 + (a**2 / (x**2 + y**2)))
    yn = y * (1 - (a**2 / (x**2 + y**2)))
    return x, y ,xn,yn

def cp (theta, x, y, epsilon, delta, alpha, a):
    C2 = (np.sin(theta) -alpha*np.cos(theta) + (delta/a + (1 + epsilon/a)*alpha) / (np.sqrt( (delta/a)**2 + (1 + epsilon/a)**2) ) )**2
    C3 = ( (x**2 + y**2)**2 ) / ((x**2 + y**2)**2 - 2* (x**2 - y**2)  + 1 )
    cp = 1 - 4*C2*C3
    return cp


R = 1
delta = 0.036
epsilon = 0.096
alpha = 2 * (np.pi / 180)
a = np.sqrt(R**2 - delta**2) - epsilon

theta = np.linspace(0, 360, 360) * (2 * np.pi / 360)
x = -epsilon + R * np.cos(theta)
y = delta + R * np.sin(theta)


plt.figure()
plt.plot(Yukowski(a, x, y)[0], Yukowski(a, x, y)[1], 'black')
plt.plot(Yukowski(a, x, y)[2], Yukowski(a, x, y)[3], 'red')
plt.axhline(0, color='gray', linewidth=1)  
plt.axvline(0, color='gray', linewidth=1)  
plt.grid(True, which='both', linestyle='--', linewidth=0.5)

plt.figure()
plt.plot(Yukowski(a, x, y)[2], cp (theta, x, y, epsilon, delta, alpha, a), 'red')
plt.gca().invert_yaxis()
plt.axhline(0, color='gray', linewidth=1)  
plt.axvline(0, color='gray', linewidth=1)  
plt.grid(True, which='both', linestyle='--', linewidth=0.5)

