# Aerodynamics
Python codes used for analysis and simulation of a steady flow around profiles and wings

The structure of this repository is given in two folders, one for the nummerical programs for profiles and another one for the wings.
Both have their instructions in order to be used:

PROFILES CODES

For a complete profile analysis, it is necessary to create a folder named "Airfoil_data" if you want to compare the numerical data calculated with the experimental/obtained data.
The files in this folder must follow the naming convention:
NACA_ProfileNumber_AngleOfAttack_ST/T (with or without trim tab).txt
After this, you must update the file name and profile variables in the code so the plots make the correct comparison.

WINGS CODES

Between lines 244 and 278, there are functions that export the calculated data to an Excel file. This is designed so that large matrices are visually clearer.
To use this code correctly, it is necessary to create two folders in the same directory as the code file: one called "Tablas" and another called "XFLR5_Data".
Tables are automatically saved with appropriate names. The XFLR5 data must be extracted from the program and assigned the name:
Wing_Grap_cd/cl/cm.txt — one file for each coefficient.
