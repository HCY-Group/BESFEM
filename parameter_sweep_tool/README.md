To generate config files, at first the required parameters that will be tested are provided in sweep.json.
Next, the code can be run in the terminal with the command python parameter_sweep.py \
    --base-config ../inputs/run_config.txt \
    --sweep sweep.json \
    --benchmark-root benchmark_runs
This will produce the folder benchmark_runs containing config files.

A slurm job can be submitted which will show the terminal outputs in the log folder and the subsequent config file used for a particular job will be at results folder.

    