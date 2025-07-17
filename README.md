# ✈️ Aerodynamics

This repository contains Python code for the analysis and simulation of steady flow around airfoils and wings. The code is structured into two main folders:

- Profiles Codes

- Wings Codes

Each folder includes specific instructions for setup and usage.

## 📄 Profiles Codes
These scripts simulate flow around airfoils and optionally compare numerical results with external data.

## Setup Instructions:

1. Create a folder named Airfoil_data in the root directory.
2. Place data files using the following naming convention:
    NACA_<ProfileNumber><AngleOfAttack>ST.txt
    NACA<ProfileNumber><AngleOfAttack>_T.txt

    - ST = no trim tab

    - T = with trim tab

### ▶️ Usage

- Open the script and modify the `filename` and `profile` variables to match the airfoil and angle of attack you're analyzing.
- Run the script to generate plots and, if applicable, compare with external data in `Airfoil_data`.

---

## 🧮 Wings Codes
Scripts for lifting-line theory and other 3D wing calculations, with optional export to Excel and comparison with XFLR5 outputs.

## Setup Instructions:

1. Create two folders in the root directory:
   - `Tablas` — for storing generated Excel tables
   - `XFLR5_Data` — for importing coefficient data exported from XFLR5

2. Place the XFLR5 output files in `XFLR5_Data` using these names:
   Wing_Grap_cd.txt
   Wing_Grap_cl.txt
   Wing_Grap_cm.txt

### 📤 Data Export

- The script automatically saves results in Excel format (in `Tablas`) for better visualization of large matrices.
- Export-related functions are located between **lines 244 and 278** in the code.
