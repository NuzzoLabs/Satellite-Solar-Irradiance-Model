# Satellite-Solar-Irradiance-Model

MATLAB orbital analysis tool developed as part of an undergraduate aerospace engineering project in 2019 for estimating solar flux incident on a satellite surface using orbital position data exported from AGI Systems Tool Kit (STK). The model processes Earth-Centered, Earth-Fixed (ECEF) position data to determine spacecraft ground track, solar incidence angles, and incident solar flux on a selected satellite face under a nadir-pointing attitude assumption.

## Features
- Reads STK Fixed Position Velocity report data
- Converts ECEF position vectors to approximate latitude and longitude
- Plots satellite ground track
- Estimates Earth-Sun distance and solar flux at Earth orbit
- Computes satellite hour angle, solar altitude, azimuth, and incidence angle
- Calculates incident solar flux on a selected spacecraft surface
- Generates diagnostic plots for solar flux, ground track, latitude, longitude, and incidence angle

## Input Data
The included sample file is:
```text
data/orbit_data_no_maneuver.tsv
```
The file is assumed to be an STK Fixed Position Velocity report with columns:
```text
Date | Time | ECEF_X_km | ECEF_Y_km | ECEF_Z_km | ECEF_VX_km_s | ECEF_VY_km_s | ECEF_VZ_km_s
```
The current model uses ECEF position columns only. Velocity columns are retained for traceability and future attitude modeling.

## Example Run Inputs
For the included dataset, use:
```text
Start of simulation day number: 170.67
End of simulation day number:   170.78
Time step:                      0.00073 days
Face:                           0
Tilt angle:                     0 degrees
```
Day number can include a decimal fraction to account for UTC time of day.

## How to Run
Open MATLAB from the repository root and run:
```matlab
addpath('src')
run_solar_flux_analysis
```
The script will prompt for the simulation start day, end day, time step, face, and tilt angle.
To run the included sample case without prompts:
```matlab
addpath('src')
results = solar_flux_analysis('data/orbit_data_no_maneuver.tsv', 170.67, 170.78, 0.00073, 0, 0, true);
```
## Repository Structure
```text
satellite-solar-flux-analysis/
├── src/
│   ├── solar_flux_analysis.m
│   └── run_solar_flux_analysis.m
├── data/
│   ├── orbit_data_no_maneuver.csv
│   └── README.md
├── docs/
│   └── original_script_2019.m
├── results/
├── README.md
└── .gitignore
```
## Methodology Notes
This model was originally developed in 2019 as an undergraduate aerospace engineering project using STK orbit data and MATLAB. It assumes a spherical Earth for latitude conversion and a nadir-pointing spacecraft attitude. The original code was refactored into a function-based MATLAB project for easier reuse and GitHub publication.

## Status
Research and educational code under active refactoring. Results should be independently validated before use in spacecraft design, mission analysis, or operational decision-making.
