# Automation of Michaelis-Menten Kinetic Analysis for Malate Dehydrogenase Datasets

Author: Emily Rose Friedman

---

## Abstract

My St. Mary’s Project explores phosphorylation as a regulatory strategy for malate dehydrogenase (MDH) via the kinetic characterization of a novel phosphomimetic mutant; S227D. Towards this enzymatic analysis I collected kinetic rate datasets for both wildtype (WT) and S227D MDH for direct comparison. These data encompass the reaction curves of the enzymes with varying concentrations of the substrate OAA, measuring the consumption of the cofactor NADH via a decrease in absorbance (340nm) over time. The custom python scripts outlined here proved successful in automating the analysis of this data according to the Michaelis-Menten kinetics model to extract relevant kinetic parameters; a process that is commonly computed by hand. This included the identification of linear regions in assay data, the fitting of linear regressions to extract V0 for each OAA concentration, conversion to specific activity values, and nonlinear curve fitting to obtain Vmax and Km values, as well as the generation of Michaelis-Menten plots for both enzymes.

---

## Introduction

For my St. Mary’s project, I have collected series of kinetic rate data for both wildtype malate dehydrogenase (MDH) and a novel phosphomimetic MDH mutant: S227D. This data encompasses the reaction curves of the enzymes with varying concentrations of substrate OAA by measuring the consumption of the cofactor NADH by decreasing absorbance (340nm) over the course of the reaction. The most common approach to statistically determine and compare the kinetic parameters of enzymes is through Michaelis-Menten Kinetics, which enzymatic reactions involving a single substrate can be assumed to follow.

The curve expressing the relationship between the initial substrate concentration [S], and the initial reaction velocity V0 has a characteristic rectangular hyperbolic shape for such enzymes, which can be stated mathematically by the Michaelis-Menten equation:

![](./michal.png)

This equation demonstrates that the initial velocity V0, the maximum velocity Vmax, and the initial substrate concentration [S] are all related quantitatively through the Michaelis constant Km. These parameters can be readily measured experimentally by running a series of enzyme assays at varying substrate concentrations and extracting the initial reaction rates (V0); which can be found in the linear portion of the data at the beginning of an assay where [S] can be regarded as constant. The parameters Vmax (the maximum rate at which an enzyme can catalyze product formation), and Km (The substrate concentration at which ½ Vmax is achieved), can then be directly derived from Michaelis-Menten plot via nonlinear regression. Nonlinear curve-fitting is the preferred method for calculation of kinetic parameters, as it avoids error propagation in lower substrate concentrations that is prevalent in a traditional double-reciprocal transformations of the data (Lineweaver-Burk plots) which are commonly employed in manual calculations of kinetic parameters.

![](./mmkinetics_schematic.png)

Manual performance of this pipeline can introduce other substantial errors. Traditionally, each assay is visually assessed for perceived linear regions; for which the range is then selected by hand. Many select the first linear region at the beginning of the assay, which is often steeper than the true initial reaction velocity. Additionally, these manual region selections do not optimize for goodness of fit when performing linear regressions, which could hurt downstream analyses such as the goodness of fit for the Michaelis-Menten curve itself. With these shortcomings in mind, I sought to optimize the Michaelis-Menten analysis pipeline by identifying low-curvature linear regions based on second derivatives to maximize goodness-of-fit, and to apply nonlinear curve-fitting to these more accurate reaction rates.

---

## Methods & Results

This pipeline consists of two primary Python programs:

- `michaelis_menten_kinetics.py`: formats and processes pathlength-corrected absorbance vs. time data, identifies linear regions of data, computes initial velocities, and converts them to specific activity which are then averaged across replicates.
- `fit_michaelis_menten.py`: Fits data to Michaelis-Menten equation, computes Vmax and Km values with standard error, and generates plots.

---

### `michaelis_menten_kinetics.py` 
[View michaelis_menten_kinetics.py](./wt500_michaelis_menten_kinetics.py)

#### Step 1. Loading and formatting data

[Wildtype data](./WT_pathcorrected_OAA_ph8_dil500.xlsx)

[Mutant data](./S227D_pathcorrected_OAA_ph8.xlsx)

- Each excel file contains 7 sheets, each corresponding to a substrate concentration.
- Load into dictionaries where each key is the substrate concentration (sheet name), and the value is a the associated dataframe:

```python
# Read excel files into dictionaries
dict_s227d = pd.read_excel(xls_s227d, sheet_name = xls_s227d.sheet_names[0:7])
dict_wt = pd.read_excel(xls_wt, sheet_name = xls_wt.sheet_names[0:7])
```

- Retain only time and absorbance data for the 3 replicates
- Rename columns to make downstream coding easier

```python
for sheet_name, df in dict_s227d.items():
    df = df.iloc[:, 1:5]  # Keep only columns 2 to 5
    df.columns = ["Time", "s227d_R1", "s227d_R2", "s227d_R3"]  # Rename columns
    dict_s227d[sheet_name] = df  # Update the dictionary

for sheet_name, df in dict_wt.items():
    df = df.iloc[:, 1:5]
    df.columns = ["Time", "wt_R1", "wt_R2", "wt_R3"]
    dict_wt[sheet_name] = df
```
---
#### Step 2. Calculate Derivatives
`def process_sheets(xls, sample_prefix):`
This function:
1.	Loops over each sheet in each Excel file.
2.	Calculates for each replicate:
- First derivative (dA/dt): Represents reaction rate over time.
- Second derivative (d²A/dt²): Measures curvature → helps find regions that are linear.

`np.gradient` numerically estimates derivatives using finite differences.

```python
# Function to compute derivatives
	# second derivative will be used to help find best linear region
def process_sheets(xls, sample_prefix):
    dict_data = pd.read_excel(xls, sheet_name=None)  # Read all sheets
    processed_dict = {}

    for sheet_name, df in dict_data.items():
        df = df.iloc[:, 1:5]  # Keep only relevant columns
        df.columns = ["Time", f"{sample_prefix}_R1", f"{sample_prefix}_R2", f"{sample_prefix}_R3"]

        # first derivative
        df[f"dA_dt_R1"] = np.gradient(df[f"{sample_prefix}_R1"], df["Time"])
        df[f"dA_dt_R2"] = np.gradient(df[f"{sample_prefix}_R2"], df["Time"])
        df[f"dA_dt_R3"] = np.gradient(df[f"{sample_prefix}_R3"], df["Time"])

        # second derivative
        df[f"d2A_dt2_R1"] = np.gradient(df[f"dA_dt_R1"], df["Time"])
        df[f"d2A_dt2_R2"] = np.gradient(df[f"dA_dt_R2"], df["Time"])
        df[f"d2A_dt2_R3"] = np.gradient(df[f"dA_dt_R3"], df["Time"])

        processed_dict[sheet_name] = df  # Store processed dataframe

    return processed_dict

# Processed data with derivatives
dict_s227d_processed = process_sheets(xls_s227d, "s227d")
dict_wt_processed = process_sheets(xls_wt, "wt")
```
---
#### Step 3. Define constants for specific activity calculations
- Used to calculate concentration from absorbance (Beer's law) and normalize for enzyme amounts in assay

```python
# Constants to be used in specific activity calculations
NADH_concentration = 0.0322  # concentration in umol of the absorbing species (c) from Beer's law

# mg of enzyme in 20ul
enzyme_amount_s227d = 0.0086666675  
enzyme_amount_wt = 0.00001592258262  # 1:500 dilution
```
