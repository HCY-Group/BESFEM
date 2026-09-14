"""
Create BESFEM config files for a multi-parameter sweep.

Example:

python parameter_sweep.py \
    --base-config run_config.txt \
    --sweep sweep.json \
    --benchmark-root benchmark_runs
"""

import argparse
import csv
import json
from datetime import datetime
from pathlib import Path


def get_arguments():
    parser = argparse.ArgumentParser(
        description="Generate BESFEM config files for a parameter sweep."
    )

    parser.add_argument(
        "--base-config",
        required=True,
        help="Path to the original BESFEM config file.",
    )

    parser.add_argument(
        "--sweep",
        required=True,
        help="Path to the sweep JSON file.",
    )

    parser.add_argument(
        "--benchmark-root",
        default="benchmark_runs",
        help="Folder where benchmark runs will be created.",
    )

    return parser.parse_args()


def read_sweep_file(sweep_path):
    with open(sweep_path, "r") as file:
        sweep = json.load(file)

    if "parameters" not in sweep:
        raise ValueError(
            "sweep.json must contain a 'parameters' object."
        )

    parameters = sweep["parameters"]

    if not isinstance(parameters, dict):
        raise ValueError(
            "'parameters' must contain parameter names and value lists."
        )

    if not parameters:
        raise ValueError(
            "The parameters object cannot be empty."
        )

    for parameter_name, parameter_values in parameters.items():
        if not isinstance(parameter_values, list):
            raise ValueError(
                f"Values for '{parameter_name}' must be a list."
            )

        if not parameter_values:
            raise ValueError(
                f"The values list for '{parameter_name}' cannot be empty."
            )

    return sweep


def parameter_exists(config_lines, parameter_name):
    for line in config_lines:
        stripped_line = line.strip()

        if not stripped_line:
            continue

        if stripped_line.startswith("#"):
            continue

        if "=" not in line:
            continue

        key, value = line.split("=", 1)

        if key.strip() == parameter_name:
            return True

    return False


def create_new_config(config_lines, selected_values):
    new_lines = []

    for line in config_lines:
        stripped_line = line.strip()

        # Keep blank lines and comments unchanged.
        if not stripped_line or stripped_line.startswith("#"):
            new_lines.append(line)
            continue

        # Keep lines without "=" unchanged.
        if "=" not in line:
            new_lines.append(line)
            continue

        key, old_value = line.split("=", 1)
        parameter_name = key.strip()

        if parameter_name in selected_values:
            new_value = selected_values[parameter_name]
            new_lines.append(
                f"{parameter_name} = {new_value}\n"
            )
        else:
            new_lines.append(line)

    return new_lines


def generate_combinations(
    parameter_names,
    parameters,
    index=0,
    current_combination=None,
):
    

    if current_combination is None:
        current_combination = {}

    if index == len(parameter_names):
        return [current_combination.copy()]

    parameter_name = parameter_names[index]
    parameter_values = parameters[parameter_name]

    combinations = []

    for parameter_value in parameter_values:
        current_combination[parameter_name] = parameter_value

        remaining_combinations = generate_combinations(
            parameter_names,
            parameters,
            index + 1,
            current_combination,
        )

        combinations.extend(remaining_combinations)

    return combinations


def main():
    args = get_arguments()

    base_config_path = Path(args.base_config)
    sweep_path = Path(args.sweep)
    benchmark_root = Path(args.benchmark_root)

    if not base_config_path.exists():
        raise FileNotFoundError(
            f"Base config was not found: {base_config_path}"
        )

    if not sweep_path.exists():
        raise FileNotFoundError(
            f"Sweep file was not found: {sweep_path}"
        )

    # Read the original BESFEM config.
    with open(base_config_path, "r") as file:
        config_lines = file.readlines()

    # Read the parameters and values to sweep.
    sweep = read_sweep_file(sweep_path)

    parameters = sweep["parameters"]
    parameter_names = list(parameters.keys())
    sweep_name = sweep.get("sweep_name", "besfem_sweep")

    for parameter_name in parameter_names:
        if not parameter_exists(config_lines, parameter_name):
            raise ValueError(
                f"Parameter '{parameter_name}' was not found "
                f"in {base_config_path}."
            )

    combinations = generate_combinations(
        parameter_names,
        parameters,
    )

    # Create a unique benchmark folder.
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S_%f")

    benchmark_directory = (
        benchmark_root / f"{timestamp}_{sweep_name}"
    )

    input_directory = benchmark_directory / "input"
    results_directory = benchmark_directory / "results"
    logs_directory = benchmark_directory / "logs"

    benchmark_directory.mkdir(
        parents=True,
        exist_ok=False,
    )

    input_directory.mkdir()
    results_directory.mkdir()
    logs_directory.mkdir()

    manifest_rows = []

    # Create one config for each parameter combination.
    for run_number, selected_values in enumerate(
        combinations,
        start=1,
    ):
        run_id = f"run_{run_number:04d}"

        config_filename = f"{run_id}.cfg"
        config_path = input_directory / config_filename

        new_config_lines = create_new_config(
            config_lines,
            selected_values,
        )

        with open(config_path, "w") as file:
            file.writelines(new_config_lines)

        manifest_row = {
            "run_id": run_id,
            "config_file": f"input/{config_filename}",
            "log_file": f"logs/{run_id}.log",
            "status": "generated",
        }

        # Add each parameter as a separate manifest column.
        for parameter_name in parameter_names:
            manifest_row[parameter_name] = selected_values[
                parameter_name
            ]

        manifest_rows.append(manifest_row)

        parameter_description = ", ".join(
            f"{name} = {value}"
            for name, value in selected_values.items()
        )

        print(
            f"Created {config_filename}: "
            f"{parameter_description}"
        )

    # Create manifest.csv.
    manifest_path = benchmark_directory / "manifest.csv"

    field_names = (
        ["run_id"]
        + parameter_names
        + [
            "config_file",
            "log_file",
            "status",
        ]
    )

    with open(manifest_path, "w", newline="") as file:
        writer = csv.DictWriter(
            file,
            fieldnames=field_names,
        )

        writer.writeheader()
        writer.writerows(manifest_rows)

    print()
    print(f"Created {len(combinations)} configurations.")
    print(f"Benchmark folder: {benchmark_directory}")
    print(f"Manifest: {manifest_path}")


if __name__ == "__main__":
    main()