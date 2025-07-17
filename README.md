# ✈️ Aerodynamics

This repository contains Python code for the analysis and simulation of steady flow around airfoils and wings. The code is structured into two main folders:

- Profiles Codes

- Wings Codes

Each folder includes specific instructions for setup and usage.

📄 Profiles Codes
These scripts simulate flow around airfoils and optionally compare numerical results with external data.

Setup Instructions:

Create a folder named Airfoil_data in the root directory.

To enable comparison with external data, add text files using the following naming convention:
NACA_<ProfileNumber>_<AngleOfAttack>_ST.txt or ..._T.txt

ST = no trim tab

T = with trim tab

Usage:

Update the file name and airfoil variables in the script to ensure accurate plotting and comparisons.

🧮 Wings Codes
These scripts analyze finite wings, exporting results and allowing for comparison with XFLR5 data.

Setup Instructions:

Create two folders in the same directory as the code:

Tablas: used to store the automatically generated Excel tables.

XFLR5_Data: place the exported data from XFLR5 here.

Data Format:

XFLR5 files must be named as follows (one file per coefficient):

Wing_Grap_cd.txt

Wing_Grap_cl.txt

Wing_Grap_cm.txt

Notes:

Data export functions are located between lines 244–278 in the script. These are designed for better readability of large matrices in Excel format.
