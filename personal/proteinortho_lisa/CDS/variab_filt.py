import pandas as pd

# Load data
data = pd.read_csv("identity_variability_cds.tsv", sep='\t')

# Calculate statistics for Variability(%)
mean_var = data['Variability(%)'].mean()
median_var = data['Variability(%)'].median()
min_var = data['Variability(%)'].min()
max_var = data['Variability(%)'].max()
std_var = data['Variability(%)'].std()

# Print statistics
print(f"Mean Variability: {mean_var:.7f}%")
print(f"Median Variability: {median_var:.7f}%")
print(f"Min Variability: {min_var:.7f}%")
print(f"Max Variability: {max_var:.7f}%")
print(f"Std Dev of Variability: {std_var:.7f}%")

# Filter threshold (Mean + 3 × Std Dev)
threshold = mean_var + 3 * std_var

# Filtering (selecting highly variable groups)
filtered_data = data[data['Variability(%)'] > threshold]
print(f"Highly variable orthogroups selected: {len(filtered_data)}")

# Save results
filtered_data.to_csv("variable_orthogroups_filt.tsv", sep='\t', index=False)
print(f"Filtered groups saved to 'variable_orthogroups_filt.tsv' (threshold: {threshold:.7f}%)")

