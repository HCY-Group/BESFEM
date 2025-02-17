import os
import numpy as np
import glob

# Define directories
no_skip_dir = "Results_3x90_0205_3"
skip_dir = "Results_3x90_0211_20"

# File types to compare
file_types = ["CnE", "CnP", "phP", "phE"]

# Number of MPI nodes
num_nodes = 6

def read_numerical_values(file_path):
    """Extract numerical values from a given text-based file."""
    values = []
    try:
        with open(file_path, "r") as file:
            for line in file:
                try:
                    values.append(float(line.strip()))
                except ValueError:
                    continue  # Ignore non-numeric lines
    except FileNotFoundError:
        print(f"Warning: File {file_path} not found!")
    return np.array(values)

# Store errors
error_summary = []

for file_type in file_types:
    for node_id in range(num_nodes):
        # Define file names
        no_skip_file = os.path.join(no_skip_dir, f"{file_type}.{'%06d' % node_id}")
        skip_file = os.path.join(skip_dir, f"{file_type}.{'%06d' % node_id}")

        # Read numerical values
        no_skip_values = read_numerical_values(no_skip_file)
        skip_values = read_numerical_values(skip_file)

        # Ensure they are the same length
        if len(no_skip_values) != len(skip_values):
            print(f"Mismatch in file sizes for {file_type} on node {node_id}")
            continue

        # Compute error metrics
        abs_error = np.abs(no_skip_values - skip_values)
        rel_error = abs_error / (np.abs(no_skip_values) + 1e-10)  # Avoid divide by zero
        l2_error = np.sqrt(np.sum(abs_error**2)) / np.sqrt(np.sum(no_skip_values**2))

        # Store results
        error_summary.append({
            "file_type": file_type,
            "node": node_id,
            "abs_error_mean": np.mean(abs_error),
            "rel_error_mean": np.mean(rel_error),
            "L2_error": l2_error
        })

# Display summary
import pandas as pd
import matplotlib.pyplot as plt

error_df = pd.DataFrame(error_summary)
print(error_df)

# Save results
error_df.to_csv("error_comparison.csv", index=False)

# Load the error data into a DataFrame
data = pd.read_csv("error_comparison.csv")

# Plot L2 error for each file type
plt.figure(figsize=(10, 5))
for file_type in data["file_type"].unique():
    subset = data[data["file_type"] == file_type]
    plt.plot(subset["node"], subset["L2_error"], marker="o", label=file_type)

plt.xlabel("Node")
plt.ylabel("L2 Error")
plt.title("L2 Error Across Nodes for Different File Types")
plt.legend()
plt.grid()
plt.show()
plt.savefig("L2_error_plot.png", dpi=300)