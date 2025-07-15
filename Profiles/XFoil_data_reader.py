# -*- coding: utf-8 -*-
"""
Created on Thu Apr 10 01:11:41 2025

@author: mdeto
"""

import numpy as np
import matplotlib.pyplot as plt

def read_airfoil_data(x):
    data = np.loadtxt(x, skiprows=3)

    x = data[:, 0]
    y = data[:, 1]
    cp = data[:, 2]
    plt.figure()
    plt.scatter(x, cp)
    plt.gca().invert_yaxis()
    plt.show()
    
    return x, y, cp

read_airfoil_data("Airfoil_data/NACA_2319_0.56_ST.txt")