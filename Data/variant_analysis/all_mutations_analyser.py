import pandas as pd

# Load the CSV file into a DataFrame
df = pd.read_csv('Data/variant_analysis/All_Variant_Frequencies/NZCP073017_allmutations.csv')

# Drop rows with NaN values in 'Variant Frequency' column
df = df.dropna(subset=['Variant Frequency'])

# Ensure 'Variant Frequency' column is numeric
df['Variant Frequency'] = pd.to_numeric(df['Variant Frequency'], errors='coerce')

# Drop rows where 'Variant Frequency' couldn't be converted to a number
df = df.dropna(subset=['Variant Frequency'])

# Remove duplicate rows
df = df.drop_duplicates()

# Check for any remaining NaNs
print(df.isna().sum())

