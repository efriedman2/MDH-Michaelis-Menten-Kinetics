#!/usr/bin/python3

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Define file path
file_dir = "~/Documents/comp_bio/3-22_OAAph8_data/reformatted_data/"

# Load excel file with average specific activity
df = pd.read_excel(f"{file_dir}avg_specific_activity_ex0.3.xlsx")

# Separate data
OAA_conc = df.iloc[:, 0].values
specific_activity_s227d = df.iloc[:, 1].values
specific_activity_wt = df.iloc[:, 2].values

###############################################################################

# Nonlinear curvefitting

# Define Michaelis-Menten equation
def michaelis_menten(S, Vmax, Km):
    return (Vmax * S) / (Km + S)

# Fit model to data
popt_s227d, _ = curve_fit(michaelis_menten, OAA_conc, specific_activity_s227d, p0=[0.08, 10])
Vmax_s227d, Km_s227d = popt_s227d

popt_wt, _ = curve_fit(michaelis_menten, OAA_conc, specific_activity_wt, p0=[25, 10])
Vmax_wt, Km_wt = popt_wt

# Print fitted equations
print(f"Fitted Michaelis-Menten equation for S227D: v = ({Vmax_s227d:.4f} * [S]) / ({Km_s227d:.4f} + [S])")
print(f"Fitted Michaelis-Menten equation for WT: v = ({Vmax_wt:.4f} * [S]) / ({Km_wt:.4f} + [S])")

# Print Vmax and Km values
print(f"S227D: Vmax = {Vmax_s227d:.2f}, Km = {Km_s227d:.2f}")
print(f"WT: Vmax = {Vmax_wt:.2f}, Km = {Km_wt:.2f}")

# Generate fine-grained points for smooth curves
S_fit = np.linspace(min(OAA_conc), max(OAA_conc), 300)
SA_fit_s227d = michaelis_menten(S_fit, Vmax_s227d, Km_s227d)
SA_fit_wt = michaelis_menten(S_fit, Vmax_wt, Km_wt)

###############################################################################

# Plot Michaelis-Menten curves

# S227D
plt.figure(figsize=(6, 5))
plt.scatter(OAA_conc, specific_activity_s227d, color='blue')
plt.plot(S_fit, SA_fit_s227d, color='blue', linestyle='--')
plt.xlabel("[OAA] (μM)")
plt.ylabel("Specific Activity (µmol NADH/min/mg)")
plt.title("Michaelis-Menten Fit (S227D)")
plt.savefig("michaelis_menten_s227d.png", dpi=300)

# WT
plt.figure(figsize=(6, 5))
plt.scatter(OAA_conc, specific_activity_wt, color='red')
plt.plot(S_fit, SA_fit_wt, color='red', linestyle='--')
plt.xlabel("[OAA] (μM)")
plt.ylabel("Specific Activity (µmol NADH/min/mg)")
plt.title("Michaelis-Menten Fit (WT)")
plt.savefig("michaelis_menten_wt.png", dpi=300)

# Combined Plot
plt.figure(figsize=(6, 5))
plt.scatter(OAA_conc, specific_activity_s227d, color='blue')
plt.scatter(OAA_conc, specific_activity_wt, color='red')
plt.plot(S_fit, SA_fit_s227d, color='blue', linestyle='--')
plt.plot(S_fit, SA_fit_wt, color='red', linestyle='--')
plt.xlabel("[OAA] (μM)")
plt.ylabel("Specific Activity (µmol NADH/min/mg)")
plt.title("Michaelis-Menten Fit (S227D & WT)")
plt.savefig("michaelis_menten_combined.png", dpi=300)
