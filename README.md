# FARGO3D–mcdust Connection

This repository provides a connection to convert gas velocity outputs from FARGO3D simulations into velocity input files compatible with the mcdust dust evolution code.

It includes:
- A standalone Python script: 'fargo3d_to_mcdust.py'
- A Jupyter notebook: 'fargo3d_to_mcdust.ipynb'
- Example simulation setups, outputs, and resulting mcdust inputs for both full and half-disk simulations

---

## Purpose

FARGO3D is a powerful hydrodynamics code used to simulate protoplanetary disks. mcdust takes gas velocity input on a cylindrical '(r, z)' grid in scale-free units. This connection performs the necessary conversion and formatting, including:

- Automatic detection of coordinate system (spherical, cylindrical, cartesian)
- Azimuthal slicing or averaging of 3D data
- Half-disk mirroring if needed
- Sound speed and scale height calculation (isothermal or adiabatic)
- Scale-free transformation
- RBF interpolation onto regular cylindrical grids
- Writing a '.inp' file usable by McDust

---

## How to Use

You can use either the Python script or the Jupyter notebook.

### 1. Edit Configuration
In 'fargo3d_to_mcdust.py' or the notebook, specify:

- 'parameter_file' and 'options_file' (optional if you have summary file in 'outputs_directory' and 'output_number' defined)
- 'outputs_directory'
- 'output_number' (or leave as 'None' to use last available)
- 'azimuthal_cut' ('None', '"average"', or integer index)
- 'mcdust_input_dir' and 'input_filename'

### 2. Run the Script
'''bash
python3 fargo3d_to_mcdust.py

Or launch the notebook and run all cells.

---

<pre>
## Repository Structure

<code>
fargo3d_mcdust_connection/
├── fargo3d_to_mcdust.py            # Main Python script
├── fargo3d_to_mcdust.ipynb         # Notebook version of the pipeline
│
├── example_fargo_simsetups/        # .par and .opt files for test simulations
│   ├── p3disof/                    # Full disk setup
│   └── p3diso/                     # Half disk setup
│
├── example_fargo_outputs/          # Gas velocity and density outputs from FARGO3D
│   ├── full_disk_sim/
│   └── half_disk_sim/
│
├── example_mcdust_inputs/          # .inp velocity files ready for McDust
│   ├── full_disk_sim.inp
│   └── half_disk_sim.inp
│
└── README.md
</code>
</pre>

---

## Requirements
Python 3.8+
numpy
scipy
matplotlib

Install using pip if needed:
'''bash
pip install numpy scipy matplotlib

---

## Example Use
To test the code, you can run the notebook with the example files provided in the example_* directories.
These include both a full-disk and a half-disk simulation, so you can see how mirroring and coordinate transformation are handled.

---

## License
None

---

## Author
Emilio Saldaña Markus:
Max Planck Institute for Solar System Research - Planetary Science Department - Planetary Formation Group
University of Texas at Austin - College of Natural Sciences - Physics Department

---

## Contributions
Feel free to use this repository by pulling it, open issues to improve the code, make suggestions to support more FARGO3D configurations or more mcdust readings, or optimize performance.
