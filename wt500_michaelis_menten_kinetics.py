#!/usr/bin/python3
import numpy as np
import pandas as pd

# Define file path
file_dir = "~/Documents/comp_bio/3-22_OAAph8_data/reformatted_data/"

# Load Excel files
xls_s227d = pd.ExcelFile(f"{file_dir}S227D_pathcorrected_OAA_ph8.xlsx")
xls_wt = pd.ExcelFile(f"{file_dir}WT_pathcorrected_OAA_ph8_dil100.xlsx")

#Read excel files into dictionaries
dict_s227d = pd.read_excel(xls_s227d, sheet_name = xls_s227d.sheet_names[0:7])
dict_wt = pd.read_excel(xls_wt, sheet_name = xls_wt.sheet_names[0:7])

#Keep only time and absorbance data, and rename columns for easier access
#Iterate through each dataframe associated with each sheetname in dictionary
for sheet_name, df in dict_s227d.items():
    df = df.iloc[:, 1:5]  # Keep only columns 2 to 5
    df.columns = ["Time", "s227d_R1", "s227d_R2", "s227d_R3"]  # Rename columns
    dict_s227d[sheet_name] = df  # Update the dictionary

for sheet_name, df in dict_wt.items():
    df = df.iloc[:, 1:5]
    df.columns = ["Time", "wt_R1", "wt_R2", "wt_R3"]
    dict_wt[sheet_name] = df

###############################################################################

# Function to clean data and compute derivatives
def process_sheets(xls, sample_prefix):
    dict_data = pd.read_excel(xls, sheet_name=None)  # Read all sheets
    processed_dict = {}

    for sheet_name, df in dict_data.items():
        df = df.iloc[:, 1:5]  # Keep only relevant columns
        df.columns = ["Time", f"{sample_prefix}_R1", f"{sample_prefix}_R2", f"{sample_prefix}_R3"]

        # Compute first derivative
        df[f"dA_dt_R1"] = np.gradient(df[f"{sample_prefix}_R1"], df["Time"])
        df[f"dA_dt_R2"] = np.gradient(df[f"{sample_prefix}_R2"], df["Time"])
        df[f"dA_dt_R3"] = np.gradient(df[f"{sample_prefix}_R3"], df["Time"])

        # Compute second derivative
        df[f"d2A_dt2_R1"] = np.gradient(df[f"dA_dt_R1"], df["Time"])
        df[f"d2A_dt2_R2"] = np.gradient(df[f"dA_dt_R2"], df["Time"])
        df[f"d2A_dt2_R3"] = np.gradient(df[f"dA_dt_R3"], df["Time"])

        processed_dict[sheet_name] = df  # Store processed dataframe

    return processed_dict

# Process both S227D and WT data
dict_s227d_processed = process_sheets(xls_s227d, "s227d")
dict_wt_processed = process_sheets(xls_wt, "wt")

################################################################################

from scipy.stats import linregress

# Constants
NADH_concentration = 0.0322  # concentration in umol of the absorbing species (c) from Beer's law

# mg of enzyme in 20ul
enzyme_amount_s227d = 0.0086666675  
enzyme_amount_wt = 0.0000796129  # 1:100 dilution

###############################################

# Function to find the best linear region for initial velocity (V₀) determination
def find_best_linear_region(df, rate_col, second_deriv_col, time_col="Time", window_size=5):
    
# Identify indices where the second derivative is among the lowest 20% of values
    threshold = np.percentile(np.abs(df[second_deriv_col]), 20)  # Dynamic threshold
    linear_region_mask = np.abs(df[second_deriv_col]) < threshold
    linear_indices = np.where(linear_region_mask)[0]

    if len(linear_indices) < window_size:
        print(f"Warning: Not enough points in the low-curvature region for {rate_col}.")
        return None, None, None, None

    # Initialize variables to track the best linear region
    best_slope, best_r2, best_start_idx, best_end_idx = None, 0, None, None

    # Iterate through possible regions
    for start_idx in range(len(linear_indices) - window_size):
        end_idx = start_idx + window_size

        time_linear = df[time_col].iloc[linear_indices[start_idx:end_idx]]
        absorbance_linear = df[rate_col].iloc[linear_indices[start_idx:end_idx]]

        # Perform linear regression
        slope, intercept, r_value, p_value, std_err = linregress(time_linear, absorbance_linear)

        # Ensure negative slope and maximize R² value
        if slope < 0 and r_value**2 > best_r2:
            best_slope, best_r2 = slope, r_value**2
            best_start_idx, best_end_idx = linear_indices[start_idx], linear_indices[end_idx - 1]

    if best_slope is not None:
        return abs(best_slope), best_start_idx, best_end_idx, best_r2
    else:
        return None, None, None, None  # No valid region found

##############################################

# Dictionary to store the average specific activity for each substrate concentration, separated by enzyme type
averaged_specific_activity_s227d = {}
averaged_specific_activity_wt = {}

# Loop through each enzyme (S227D and WT)
for enzyme_dict, enzyme_name, enzyme_amount, averaged_specific_activity in zip(
    [dict_s227d_processed, dict_wt_processed], 
    ["S227D", "WT"], 
    [enzyme_amount_s227d, enzyme_amount_wt], 
    [averaged_specific_activity_s227d, averaged_specific_activity_wt]
):
    for sheet_name, df in enzyme_dict.items():
        substrate_concentration = sheet_name  # Assuming sheet names represent substrate concentrations
        specific_activities = []
        
        # Loop through each replicate (R1, R2, R3)
        for replicate in [f"{enzyme_name.lower()}_R1", f"{enzyme_name.lower()}_R2", f"{enzyme_name.lower()}_R3"]:
            # Ensure that the second derivative column exists
            second_deriv_col = f"d2A_dt2_{replicate[-2:]}"
            
            if second_deriv_col not in df.columns:
                print(f"Error: Column {second_deriv_col} not found in {sheet_name} for {replicate}")
                continue
            
            # Step 1: Find the best linear region for V₀ determination
            best_slope, best_start_idx, best_end_idx, min_second_derivative = find_best_linear_region(
                df,
                replicate,  # First derivative column (e.g., "s227d_R1", "wt_R1", etc.)
                f"d2A_dt2_{replicate[-2:]}"  # Second derivative column (e.g., "d2A_dt2_R1", "d2A_dt2_R2", etc.)
            )
            
            # Step 2: Calculate V₀ in abs/sec (slope from linear region)
            if best_slope is not None:
                V0 = abs(best_slope)  # The magnitude of the slope gives the initial velocity (abs/sec)
                print(f"V₀ for {replicate} in {sheet_name}: {V0:.6f} abs/sec")
                print(f"Linear region: Time from {df['Time'].iloc[best_start_idx]} to {df['Time'].iloc[best_end_idx]} (seconds)")
                
                # Step 3: Convert V₀ to specific activity (µmol NADH/min/mg)
                enzyme_units = ((V0 * 60) * NADH_concentration) # convert abs/sec to abs/min before multiplying by c
                specific_activity = (enzyme_units / enzyme_amount)
                
                print(f"Specific Activity for {replicate} in {sheet_name}: {specific_activity:.4f} µmol NADH/min/mg")
                
                # Step 4: Store specific activity for averaging
                specific_activities.append(specific_activity)
            else:
                print(f"No valid linear region found for {replicate} in {sheet_name}")
        
        # Step 5: Calculate and store the average specific activity for this substrate concentration
        if specific_activities:
            averaged_specific_activity[substrate_concentration] = np.mean(specific_activities)
        else:
            averaged_specific_activity[substrate_concentration] = None

# Create a DataFrame for the averaged specific activities
df_s227d = pd.DataFrame(list(averaged_specific_activity_s227d.items()), columns=["Substrate Concentration", "Average Specific Activity (S227D)"])
df_wt = pd.DataFrame(list(averaged_specific_activity_wt.items()), columns=["Substrate Concentration", "Average Specific Activity (WT)"])

# Merge the S227D and WT DataFrames on the "Substrate Concentration" column
df_combined = pd.merge(df_s227d, df_wt, on="Substrate Concentration", how="outer")

# Save the DataFrame to a CSV or Excel file
df_combined.to_csv("averaged_specific_activity.csv", index=False)
# or to save as Excel
# df_combined.to_excel("averaged_specific_activity.xlsx", index=False)

# Print the result
print("Averaged Specific Activities saved successfully!")
