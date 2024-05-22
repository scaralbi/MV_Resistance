import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from matplotlib.colors import to_rgba

# Define the light cycle pattern (dark and light cycles)
light_cycle_pattern = [0]*40 + [50]*10 + [0]*10 + [80]*10 + [0]*10 + [120]*10 + [0]*10 + \
                      [170]*10 + [0]*10 + [220]*10 + [0]*10 + [520]*10 + [0]*10 + [900]*10 + [0]*10 + [1200]*10 + [0]*20

# Load the 'All' sheet from the Excel file
path_to_excel_file = '/Volumes/Albi/Cambridge/PhD/Data/Echem_MV/Oxygen_Electrodes/Oxygen_Electrodes_Data_All.xlsx'
all_data_sheet = pd.read_excel(path_to_excel_file, sheet_name='All')
all_data_sheet = all_data_sheet.iloc[1:] # Remove the first row containing redundant labels

# Reshape the data and calculate the average and standard deviation
reshaped_data = []
for i in range(0, all_data_sheet.shape[1], 4):
    condition = all_data_sheet.columns[i].strip('"')
    time_col = all_data_sheet.columns[i]
    for j in range(1, 4):
        replicate_col = all_data_sheet.columns[i + j]
        replicate_data = all_data_sheet[[time_col, replicate_col]].copy()
        replicate_data.columns = ["Time", "Voltage"]
        replicate_data["Voltage"] = pd.to_numeric(replicate_data["Voltage"], errors='coerce').fillna(0) # Handle non-numeric values
        replicate_data["Condition"] = condition
        replicate_data["Replicate"] = j
        reshaped_data.append(replicate_data)

reshaped_data_df = pd.concat(reshaped_data, ignore_index=True)
average_std_dev_data = reshaped_data_df.groupby(["Condition", "Time"]).agg(
    Average_Voltage=("Voltage", "mean"),
    Standard_Deviation=("Voltage", "std")
).reset_index()
average_std_dev_data["Time"] = pd.to_numeric(average_std_dev_data["Time"])

# Example x-offsets for manual alignment (ADJUST THESE VALUES AS NEEDED)
example_x_offsets = {
    "WTH-MV": 0,
    "WTH+MV": 0,
    "WTN+MV": 0,
    "MV6-MV": 0,
    "MV6+MV": 0,
}

def visualize_voltage(data, x_offsets, light_cycle_pattern, max_time, title="Voltage vs Time"):
    # Split data based on condition
    mv_neg_data = data[data["Condition"].str.contains("-MV")]
    mv_pos_data = data[data["Condition"].str.contains(r"\+MV")]

    fig, axes = plt.subplots(1, 2, figsize=(15, 7))

    for ax, condition_data, title in zip(axes, [mv_neg_data, mv_pos_data], ["-MV Conditions", "+MV Conditions"]):
        # Plot the voltage data for each condition with the specified x-offsets
        for condition in condition_data["Condition"].unique():
            condition_specific_data = condition_data[condition_data["Condition"] == condition]
            time_values = condition_specific_data["Time"] + x_offsets.get(condition, 0)
            average_voltage = condition_specific_data["Average_Voltage"]
            std_dev = condition_specific_data["Standard_Deviation"]
            ax.plot(time_values, average_voltage, label=condition)
            ax.fill_between(time_values, average_voltage - std_dev, average_voltage + std_dev, alpha=0.2)
        for i in range(len(light_cycle_pattern)):
            color_intensity = light_cycle_pattern[i] / max(light_cycle_pattern)
            color = to_rgba('black', 1 - color_intensity) # Black color with varying alpha value based on intensity
            if i % 10 == 0:
                ax.axvspan(i, i + 10, facecolor=color, alpha=0.2)
        # Add a legend and labels
        ax.legend(loc='upper left')
        ax.set_title(title)
        ax.set_xlabel("Time (Minutes)")
        ax.set_ylabel("Voltage")
        ax.set_xlim(0, max_time) # Set the maximum time limit
    plt.show()

# Visualize the data with the example x-offsets
visualize_voltage(average_std_dev_data, example_x_offsets, light_cycle_pattern, max_time=250, title="Voltage vs Time")



# Constants for converting voltage to oxygen concentration
max_oxygen = 236 # Maximum oxygen in nmol O2 mL-1
delta_voltage = max_oxygen # Assuming max steady state voltage corresponds to max_oxygen and electric voltage = 0
conversion_factor = max_oxygen / delta_voltage # Conversion factor in nmol O2 mL-1 mV-1

# Modify the reshaped_data to include the oxygen concentration instead of voltage
for replicate_data in reshaped_data:
    replicate_data["Concentration"] = replicate_data["Voltage"] * conversion_factor
    replicate_data.drop(columns=["Voltage"], inplace=True)

# Concatenate reshaped_data
reshaped_data_df = pd.concat(reshaped_data, ignore_index=True)

# Define a function to get the color for each condition
def get_color(condition, color_map=plt.cm.viridis):
    unique_conditions = sorted(list(set(reshaped_data_df["Condition"])))
    condition_idx = unique_conditions.index(condition)
    return color_map(condition_idx / len(unique_conditions))

def visualize_concentration(data, x_offsets, light_cycle_pattern, max_time, title="Concentration vs Time"):
    fig, axes = plt.subplots(1, 2, figsize=(15, 7))

    for ax, mv_condition, title in zip(axes, ["-MV", r"+MV"], ["-MV Conditions", "+MV Conditions"]):
        # Plot the concentration data for each condition and replicate
            for condition in data["Condition"].unique():
                if mv_condition not in condition:
                    continue
            condition_specific_data = data[data["Condition"] == condition]
            color = 'black' # Use black color for scatter plot
            for replicate in condition_specific_data["Replicate"].unique():        
                replicate_data = condition_specific_data[condition_specific_data["Replicate"] == replicate]
                time_values = replicate_data["Time"] + x_offsets.get(condition, 0)
                concentration_values = replicate_data["Concentration"]
                ax.scatter(time_values, concentration_values, label=condition if replicate == 1 else "", color=color, s=5)
            # ... Same code as before ...
            for i in range(len(light_cycle_pattern)):
                color_intensity = light_cycle_pattern[i] / max(light_cycle_pattern)
                color = to_rgba('black', 1 - color_intensity) # Black color with varying alpha value based on intensity
                if i % 10 == 0:
                    ax.axvspan(i, i + 10, facecolor=color, alpha=0.2)

                # Add a legend and labels
                ax.legend(loc='upper left')
                ax.set_title(title)
                ax.set_xlabel("Time (Minutes)")
                ax.set_ylabel("[Oxygen] (nmol/ml)")
                ax.set_xlim(0, max_time) # Set the maximum time limit
                
    plt.show()

# Visualize the data with the example x-offsets
visualize_concentration(reshaped_data_df, example_x_offsets, light_cycle_pattern, max_time=250, title="Concentration vs Time")


# Calculating the maximum voltage using the concentration column
max_voltage = reshaped_data_df["Concentration"].max() * delta_voltage / max_oxygen

# Conversion factor now using max_voltage
conversion_factor = max_oxygen / max_voltage

# Function to detect regimes and calculate the steady state max
def detect_regimes(data, window_size=10):
    regimes = []
    for condition in data["Condition"].unique():
        condition_specific_data = data[data["Condition"] == condition]
        for replicate in condition_specific_data["Replicate"].unique():
            replicate_data = condition_specific_data[condition_specific_data["Replicate"] == replicate]
            # You've already calculated the "Concentration" earlier, so no need to recalculate it
            # Calculating the derivative to find changes in concentration
            derivative = np.gradient(replicate_data["Concentration"], replicate_data["Time"])
            # Finding peaks in the derivative to detect focal points
            peaks, _ = find_peaks(derivative, height=0)
            # Analyzing regimes around the peaks
            for peak in peaks:
                regime_data = replicate_data.iloc[peak - window_size:peak + window_size]
                steady_state_max = regime_data["Concentration"].max()
                regimes.append((condition, replicate, steady_state_max))
    return regimes


def find_light_cycle_patterns(data, window_size=10):
    light_cycle_regimes = []

    for condition in data["Condition"].unique():
        condition_specific_data = data[data["Condition"] == condition]
        for replicate in condition_specific_data["Replicate"].unique():
            replicate_data = condition_specific_data[condition_specific_data["Replicate"] == replicate]
            # Calculate the derivative to find changes in concentration
            derivative = np.gradient(replicate_data["Concentration"], replicate_data["Time"])
            # Finding peaks in the derivative to detect focal points
            peaks, _ = find_peaks(derivative, height=0)
            # Analyzing regimes around the peaks
            for peak in peaks:
                regime_data = replicate_data.iloc[peak - window_size:peak + window_size]
                light_cycle_regimes.append((condition, replicate, peak, regime_data))

    return light_cycle_regimes

light_cycle_regimes = find_light_cycle_patterns(reshaped_data_df)


def visualize_regimes(data, regimes, light_cycle_pattern, max_time, title="Regimes vs Time"):
    fig, ax = plt.subplots(figsize=(15, 7))
    # Plot original data
    for condition in data["Condition"].unique():
        condition_specific_data = data[data["Condition"] == condition]
        time_values = condition_specific_data["Time"]
        concentration_values = condition_specific_data["Concentration"]
        ax.scatter(time_values, concentration_values, label=condition, s=5)
    # Plot detected regimes
    for condition, replicate, steady_state_max in regimes:
        ax.axhline(steady_state_max, color='red', linestyle='--')
    # Plot light-dark cycles (in greyscale)
    for i in range(len(light_cycle_pattern)):
        color_intensity = light_cycle_pattern[i] / max(light_cycle_pattern)
        color = to_rgba('black', 1 - color_intensity)
        if i % 10 == 0:
            ax.axvspan(i, i + 10, facecolor=color, alpha=0.2)
    # Add labels and legend
    ax.set_title(title)
    ax.set_xlabel("Time (Minutes)")
    ax.set_ylabel("[Oxygen] (nmol/ml)")
    ax.set_xlim(0, max_time)
    ax.legend(loc='upper left')
    plt.show()

visualize_regimes(reshaped_data_df, regimes, light_cycle_pattern, max_time=250, title="Regimes vs Time")
