#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
FARGO3D Output to mcdust Input Converter
----------------------------------------

This script reads the output of a 3D FARGO3D gas simulation and converts
the velocity data into a scale-free cylindrical grid compatible with the
mcdust dust evolution code.

Features:
- Automatically reads parameters from '.par', '.opt', or summary files
- Detects coordinate system (spherical, cylindrical, cartesian)
- Converts to cylindrical (r, z) and handles half-disk mirroring
- Computes sound speed and scale height (isothermal or adiabatic)
- Applies scale-free transformation
- Performs RBF interpolation onto a regular grid
- Outputs McDust-compatible '.inp' file

Author: Emilio Saldaña Markus
Project: Connecting FARGO3D and mcdust to Study Dust Evolution in Structured
Protoplanetary Disks
Institutions: Max Planck Institute for Solar System Research - Planetary
Science Department - Planetary Formation Group
Date: August 2025
Repository: https://github.com/emiliosaldanam/fargo3d_mcdust_connection

Requirements:
- Python 3.8+
- numpy
- matplotlib
- scipy

Usage:
- Configure the parameter, options, output, and input settings in the script
- Run: python3 fargo3d_to_mcdust.py

"""

# Imports
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import RBFInterpolator

# ----- EDIT HERE: -----

# State directories and files of FARGO3D simulation and mcdust input

# Write .par and .opt files (optional). Leave as None if want to read
# parameters and options from summary file.
# NOTE: if no parameter file inputed, MUST input output_number desired.
parameter_file       = None  # Write .par file as str
options_file         = None  # Write .opt file as str

# Directory where outputs are located (including summary file).
outputs_directory    = ""  # Write directory where sim outputs are located

# Directory to put mcdust file created and name of file.
mcdust_input_dir     = ""  # Write directory of new mcdust input file
input_filename       = ""  # Write name of new mcdust input file

# Select options of data transformation. Edit if want a specific option.

output_number        = None  # Default (None) is last output
mcdust_grid_size     = None  # Default (None) is 128 (square grid)
azimuthal_cut        = None  # None is Nx/2, value (other cut), "average"

# ----- STOP EDITING HERE -----

# Considering no .opt and .par file:
if not options_file:
    summary_file = f"{outputs_directory}/summary{output_number}.dat"
    options_file = summary_file


# Read the information of the FARGO 3D simulation to consider parameters.

# Function to read parameters
class Parameters:
    """
    Class for reading simulation parameters from a parameter file.

    Parameters
    ------
    paramfile : str
        Name of the parameter file, normally "variables.par".

    Attributes
    ---------
    par : dict
        Dictionary with parameter names as keys and their corresponding values.
    """

    def __init__(self, paramfile):
        self.par = {}

        # Reading of parameter file.
        try:
            with open(paramfile, 'r') as file:
                lines = file.readlines()
        except IOError:
            print(f"Error: {paramfile} not found.")
            return

        # Extracting parameters and setting as attributes.
        for line in lines:
            line = line.strip()
            if not line or line.startswith("#"):  # Skip empty or comment lines
                continue

            tokens = line.split()
            if len(tokens) < 2:
                print(f"Warning: Skipping malformed line: '{line}'")
                continue

            name, value = tokens[0], tokens[1]
            value = self._cast_value(value)
            self.par[name] = value
            setattr(self, name.lower(), value)

    def _cast_value(self, value):
        """
        Attempts to cast value to float, int, or leave as string.
        """
        try:
            if '.' in value or 'e' in value.lower():
                return float(value)
            else:
                return int(value)
        except ValueError:
            return value  # Keeps as string

    def print_all_parameters(self):
        """
        Print all parameters stored in self.par.
        """
        print("Stored Simulation Parameters:\n")
        for key, value in self.par.items():
            print(f"{key}: {value}")

# Function to read summary file in case of no parameter file
class SummaryParameters:
    """
    Class for reading simulation parameters from a FARGO3D summary file.

    Parameters
    ----------
    summary_file : str
        Name of the summary file, typically like 'summary400.dat'.

    Attributes
    ----------
    par : dict
        Dictionary with parameter names as keys and their corresponding values.
    """

    def __init__(self, summary_file):
        self.par = {}

        try:
            with open(summary_file, 'r') as file:
                lines = file.readlines()
        except IOError:
            print(f"Error: {summary_file} not found.")
            return

        for line in lines:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("***"):
                continue

            # Try tab split first, then fallback to space
            if "\t" in line:
                tokens = line.split("\t")
            else:
                tokens = line.split()

            if len(tokens) < 2:
                continue  # silently skip malformed lines

            name, value = tokens[0], tokens[1]
            value = self._cast_value(value)
            self.par[name] = value
            setattr(self, name.lower(), value)

    def _cast_value(self, value):
        """
        Attempts to cast value to float, int, or leave as string.
        """
        try:
            if '.' in value or 'e' in value.lower():
                return float(value)
            else:
                return int(value)
        except ValueError:
            return value

    def print_all_parameters(self):
        """
        Print all parameters stored in self.par.
        """
        print("Stored Simulation Parameters:\n")
        for key, value in self.par.items():
            print(f"{key}: {value}")

if parameter_file:
    params = Parameters(parameter_file)
    if output_number is None:
        output_number = (params.ntot - 1) // params.ninterm
else:
    params = SummaryParameters(summary_file)


# Creation of grid and data loading considering type of coordinate system.

# Helper functions:
def detect_coordinate_system(opt_file_path):
    """
    Detects the coordinate system from a .opt file based on defined flags.
    """
    with open(opt_file_path, 'r') as f:
        content = f.read().upper()
        if "-DSPHERICAL" in content:
            return "spherical"
        elif "-DCYLINDRICAL" in content:
            return "cylindrical"
        elif "-DCARTESIAN" in content:
            return "cartesian"
        else:
            raise ValueError("No known coordinate system flag found in .opt file.")

def spherical_to_cylindrical(r, theta, vr, vtheta, vphi):
    """
    Converts spherical coordinates and velocities to cylindrical.
    """
    v_rho = vr * np.sin(theta) + vtheta * np.cos(theta)
    v_z = vr * np.cos(theta) - vtheta * np.sin(theta)
    v_phi = vphi
    return v_rho, v_phi, v_z

def cartesian_to_cylindrical(x, z, vx, vz, vy):
    """
    Converts Cartesian positions and velocities to cylindrical coordinates.
    Assumes:
        x: horizontal
        z: horizontal
        y: vertical
    """
    r_cyl = np.sqrt(x**2 + z**2)
    phi = np.arctan2(z, x)

    v_rho = (x * vx + z * vz) / r_cyl
    v_phi = (-z * vx + x * vz) / r_cyl
    v_z = vy
    return r_cyl, phi, v_rho, v_phi, v_z

def reduce_phi(data: np.ndarray, phi_cut = None) -> np.ndarray:
    """
    Reduce a 3D (Nz, Ny, Nx) array along the φ axis according to user
    specification.

    Parameters
    ----------
    data : np.ndarray
        The input 3D array to reduce.
    phi_cut : int, str, or None
        - If None, returns middle φ slice if Nx > 1, else φ=0.
        - If int, returns φ slice at that index.
        - If "average", returns azimuthal average over φ (axis=2).

    Returns
    -------
    np.ndarray
        2D array reduced along φ axis, shape (Nz, Ny).
    """
    Nz, Ny, Nx = data.shape

    if Nx == 1:
        return data[:, :, 0]
    elif phi_cut == "average":
        return np.mean(data, axis=2)
    elif isinstance(phi_cut, int):
        return data[:, :, phi_cut]
    else:
        return data[:, :, Nx // 2]

def load_vel_component(name, output_number=1, phi_cut = None):
    """
    Loads one gas velocity component at a specific output and applies
    azimuthal reduction.

    Parameters
    ----------
    name : str
        Name of the velocity component (e.g., "gasvy").
    output_number : int, optional
        Output number to load.
    phi_cut : int or str or None, optional
        - If None (default), cuts at middle φ index if Nx > 1, or 0 if Nx == 1.
        - If int, uses that index as the φ cut.
        - If "average", returns azimuthal average over φ.

    Returns
    -------
    2D ndarray
        Data on (z, r) grid after specified φ cut or average.
    """
    fname = f"{outputs_directory}/{name}{output_number}.dat"
    data = np.fromfile(fname).reshape((Nz, Ny, Nx))
    return reduce_phi(data, phi_cut)

# Detecting coordinate system and loading parameters and velocity components
coord_system = detect_coordinate_system(options_file)
Nx = params.nx
Ny = params.ny
Nz = params.nz
vx = load_vel_component("gasvx", output_number, azimuthal_cut)
vy = load_vel_component("gasvy", output_number, azimuthal_cut)
vz = load_vel_component("gasvz", output_number, azimuthal_cut)

# Loading and building the appropriate coordinate grid (after Nx/Ny/Nz are known)
if coord_system == "spherical":
    r_min = params.ymin
    r_max = params.ymax
    theta_min = params.zmin
    theta_max = params.zmax

    r = np.linspace(r_min, r_max, Ny)  # Creates 1D array of r values
    theta = np.linspace(theta_min, theta_max, Nz)  # Creates 1D array of theta values
    Theta, R = np.meshgrid(theta, r, indexing="ij")  # Creates 2D spherical grid (Ny, Nz) (r, theta)

    # Convert to cylindrical
    Z = R * np.cos(Theta)
    r_cyl = R * np.sin(Theta)

    v_rho, v_phi, v_z = spherical_to_cylindrical(R, Theta, vx, vz, vy)

elif coord_system == "cylindrical":
    r_min = params.ymin
    r_max = params.ymax
    z_min = params.zmin
    z_max = params.zmax

    r_cyl = np.linspace(r_min, r_max, Ny)
    Z = np.linspace(z_min, z_max, Nz)

    v_rho = vx
    v_phi = vy
    v_z = vz

elif coord_system == "cartesian":
    x_min = params.ymin
    x_max = params.ymax
    z_min = params.zmin
    z_max = params.zmax

    x = np.linspace(x_min, x_max, Ny)
    z = np.linspace(z_min, z_max, Nz)
    X, Z = np.meshgrid(x, z, indexing="ij")

    r_cyl, phi, v_rho, v_phi, v_z = cartesian_to_cylindrical(X, Z, vx, vz, vy)

    print("⚠ Cartesian to cylindrical conversion applied. Verify axis assumptions.")

else:
    raise ValueError(f"Unknown coordinate system: {coord_system}")


# Checks if Fargo 3D simulation is half or full disk. Mirrors disk of half disk.

# Function to check if full or half disk.
def is_half_disk(options_file):
    """
    Checks whether the simulation uses the -DHALFDISK option in the .opt file.

    Parameters:
    -----------
    options_file : str
        Full path to the .opt file.

    Returns:
    --------
    bool
        True if -DHALFDISK is enabled, False otherwise.
    """
    if not os.path.exists(options_file):
        raise FileNotFoundError(f".opt file not found at: {options_file}")

    with open(options_file, "r") as f:
        content = f.read()

    return "-DHALFDISK" in content

# Function to mirror the half disk if necessary.
def mirror_half_disk(r_cyl, Z, v_rho, v_z):
    """
    Mirrors the upper half-disk data to produce a full disk including the midplane.

    Assumes input arrays have shape (Nz, Nr) = (z, r),
    and mirrors vertically across the midplane including the midplane.

    Returns
    -------
    r_cyl_full : ndarray, shape (2*Nz, Nr)
    Z_full : ndarray, shape (2*Nz, Nr)
    v_rho_full : ndarray, shape (2*Nz, Nr)
    v_z_full : ndarray, shape (2*Nz, Nr)
    """
    # Flip full arrays
    Z_mirror = -Z[::-1, :]
    r_cyl_mirror = r_cyl[::-1, :]          # flipped for consistency
    v_rho_mirror = v_rho[::-1, :]          # even symmetry
    v_z_mirror = -v_z[::-1, :]             # odd symmetry

    # Concatenate mirrored + original (bottom + top)
    Z_full = np.concatenate((Z_mirror, Z), axis = 0)
    r_cyl_full = np.concatenate((r_cyl_mirror, r_cyl), axis = 0)
    v_rho_full = np.concatenate((v_rho_mirror, v_rho), axis = 0)
    v_z_full = np.concatenate((v_z_mirror, v_z), axis = 0)

    return r_cyl_full, Z_full, v_rho_full, v_z_full

# Function to mirror other fields when necessary.
def mirror_scalar_field(field):
    """
    Mirror a scalar field (e.g., cs or gas_density) across the disk midplane,
    including the midplane itself.

    Assumes input shape (Nz, Nr). Output shape will be (2*Nz, Nr).
    """
    field_mirror = field[::-1, :]  # Full vertical flip.
    return np.concatenate((field_mirror, field), axis = 0)

if is_half_disk(options_file):
    r_cyl, Z, v_rho, v_z = mirror_half_disk(r_cyl, Z, v_rho, v_z)
    print("Half disk detected, disk mirrored")


# Scale-free transformation and interpolation.

# Helper functions:
def is_isothermal(options_file: str) -> bool:
    """
    Check if the simulation uses an isothermal disk setup.

    Parameters
    ----------
    options_file : str
        Path to the options file used in the simulation.

    Returns
    -------
    bool
        True if the disk is isothermal, False otherwise.
    """
    with open(options_file, "r") as f:
        return "-DISOTHERMAL" in f.read()

def load_gas_energy(filepath: str, shape: tuple[int, int, int]) -> np.ndarray:
    """
    Load the gas energy file and reshape it.

    Parameters
    ----------
    filepath : str
        Full path to the gas energy file.
    shape : tuple
        Shape to reshape the binary data (Nz, Ny, Nx).

    Returns
    -------
    np.ndarray
        Reshaped gas energy array.
    """
    return np.fromfile(filepath).reshape(shape)

def load_density(filepath: str, shape: tuple[int, int, int]) -> np.ndarray:
    """
    Load the gas density file and reshape it.

    Parameters
    ----------
    filepath : str
        Full path to the gas density file.
    shape : tuple
        Shape to reshape the binary data (Nz, Ny, Nx).

    Returns
    -------
    np.ndarray
        Reshaped gas density array.
    """
    return np.fromfile(filepath).reshape(shape)

def compute_sound_speed(
    params,
    outputs_directory: str,
    output_number: int,
    options_file: str,
    shape: tuple[int, int, int],
    phi_cut = None
) -> np.ndarray:
    """
    Compute the sound speed (cs) based on the simulation type (isothermal or adiabatic).

    Parameters
    ----------
    params : object
        Object containing simulation parameters (must have .aspectratio and .gamma).
    outputs_directory : str
        Directory where the output files are stored.
    output_number : int
        Output number to read from.
    options_file : str
        File that contains options to determine if the disk is isothermal.
    shape : tuple
        Tuple with the shape (Nz, Ny, Nx) to reshape binary data.
    phi_cut : int or str or None
        Azimuthal cut method: "average", index, or None for default behavior.

    Returns
    -------
    np.ndarray
        2D sound speed array in (z, r) plane.
    """
    energy_path = f"{outputs_directory}/gasenergy{output_number}.dat"
    gas_energy_3d = load_gas_energy(energy_path, shape)
    gas_energy = reduce_phi(gas_energy_3d, phi_cut)

    if is_isothermal(options_file):
        return gas_energy

    gamma = params.gamma
    density_path = f"{outputs_directory}/gasdens{output_number}.dat"
    gas_density_3d = load_density(density_path, shape)
    gas_density = reduce_phi(gas_density_3d, phi_cut)

    return np.sqrt(gamma * (gamma - 1) * gas_energy / gas_density)

def compute_scale_height(cs: np.ndarray, r_cyl: np.ndarray, G: float = 1.0, M: float = 1.0) -> np.ndarray:
    """
    Compute the gas scale height H_g = c_s / Omega, where Omega is the Keplerian angular velocity.

    Parameters
    ----------
    cs : np.ndarray
        Sound speed array (same shape as r_cyl).
    r_cyl : np.ndarray
        Cylindrical radius array.
    G : float, optional
        Gravitational constant. Default is 1 (code units).
    M : float, optional
        Central mass. Default is 1 (code units).

    Returns
    -------
    np.ndarray
        Scale height array H_g.
    """
    omega = np.sqrt((G * M) / (r_cyl ** 3))
    return cs / omega

# Computing sound speed considering the type of simulation
cs = compute_sound_speed(params, outputs_directory, output_number, options_file, (Nz, Ny, Nx))

# Mirror cs if needed
if is_half_disk(options_file):
    cs = mirror_scalar_field(cs)

Hg = compute_scale_height(cs, r_cyl)

# Make data scale free
v_rho_sf = v_rho / cs
Z_sf = Z / Hg


# Interpolation of data

# Interpolation parameters:

# Number of points in regular cylindrical grid (Square grid always)
if mcdust_grid_size is None:
    mcdust_grid_size = 128
N_r_interp = mcdust_grid_size
N_z_interp = mcdust_grid_size

# Radial interpolation range in AU
r_min_interp = r_cyl.min()
r_max_interp = r_cyl.max()

# Vertical interpolation range in AU
z_min_interp = Z_sf.min()
z_max_interp = Z_sf.max()

# Interpolation:

# Creating regular cylindrical grid for interpolation (mcdust grid)
r_interp = np.linspace(r_min_interp, r_max_interp, N_r_interp)  # Creates radial 1D array for mcdust
z_interp = np.linspace(z_min_interp, z_max_interp, N_z_interp)  # Creates z 1D array for mcdust
RR_interp, ZZ_interp = np.meshgrid(r_interp, z_interp, indexing = "ij")  # Creates 2D array of R and Z for mcdust

# Preparing input points for interpolation
points = np.column_stack((r_cyl.flatten(), Z.flatten()))  # Use scale free coordinates.

# Interpolate v_rho / c_s (scale-free)
values_vrho = v_rho_sf.flatten()

rbf_vrho = RBFInterpolator(points, values_vrho, neighbors = 50, smoothing = 0.0)
v_rho_interp = rbf_vrho(np.column_stack((RR_interp.flatten(), ZZ_interp.flatten())))
v_rho_interp = v_rho_interp.reshape(RR_interp.shape)

# Interpolate v_z (unscaled, not used in mcdust yet)
values_vz = v_z.flatten()

rbf_vz = RBFInterpolator(points, values_vz, neighbors = 50, smoothing = 0.0)
v_z_interp = rbf_vz(np.column_stack((RR_interp.flatten(), ZZ_interp.flatten())))
v_z_interp = v_z_interp.reshape(RR_interp.shape)


# Writing of mcdust compatible file to run an mcdust simulation

os.makedirs(mcdust_input_dir, exist_ok = True)

filename = f"{mcdust_input_dir}/{input_filename}.inp"
with open(filename, "wb") as f:
    np.savetxt(f, [len(r_interp)], fmt="%d")
    f.write(b"#r,z(AU)\n")
    np.savetxt(f, np.c_[r_interp.flatten(), z_interp.flatten()])  # r and z_sf
    f.write(b"#RR (AU), v_rho (cs)\n")
    np.savetxt(f, np.c_[RR_interp.flatten(), v_rho_interp.flatten()])

print(f"Saved mcdust radial velocity input to {filename}")
