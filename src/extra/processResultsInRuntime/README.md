Watch and run in the `rub.pbs` script

```bash
# One-time install (no root needed)
pip install watchdog --user

# Run alongside your simulation
module load ParaView/6.0.1-foss-2025b
python watch_and_render.py ./results ./screenshots &
WATCHER_PID=$!

./my_fluid_simulation --output ./results

# After sim finishes, let the queue drain then stop
wait $WATCHER_PID
```

Example configuration file `render_config.toml`

```ini
# render_config.toml
# Place this next to your .pvsm files.
# All paths are relative to the directory containing this file,
# unless they start with / (absolute).

results_dir     = "../results"
screenshots_dir = "../screenshots"

pvpython        = "pvpython"
settle_time     = 2.0
worker_threads  = 2

# VTK filename patterns to watch for (Python regex)
patterns = [
    'fluid_[\d.]+_particles\.vtk$',
    'boundary_[\d.]+_particles\.vtk$',
]

# label = "path/to/state.pvsm"  (relative to this file's directory)
[states]
pressure = "results_pressure_vertical_legend.pvsm"
velocity = "results_velocity_vertical_legend.pvsm"
```

Example usage

```bash
# Standard use — config file does everything
python watch_and_render.py --config ./sim_A/render_config.toml

# Override results dir at runtime (e.g. different HPC scratch path)
python watch_and_render.py --config ./render_config.toml \
    --results /scratch/$SLURM_JOB_ID/results

# Completely CLI-driven, no config file
python watch_and_render.py \
    --results ./results \
    --screenshots ./shots \
    --script-dir ./states \
    --states pressure:results_pressure.pvsm velocity:results_velocity.pvsm \
    --patterns "fluid_[\d.]+_particles\.vtk$" "boundary_[\d.]+_particles\.vtk$"

# Different problem — just a different config file, same script
python watch_and_render.py --config ./sim_B/render_config.toml
```
