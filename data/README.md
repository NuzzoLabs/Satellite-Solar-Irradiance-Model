## Input Data
`orbit_data_no_maneuver.tsv` is sample STK Fixed Position Velocity report data exported without headers.
Assumed columns:
1. Date
2. Time / elapsed timestamp
3. ECEF X position, km
4. ECEF Y position, km
5. ECEF Z position, km
6. ECEF X velocity, km/s
7. ECEF Y velocity, km/s
8. ECEF Z velocity, km/s
   
The MATLAB model currently uses only columns 3-5 for ECEF position. Velocity columns are retained because they came from the STK Fixed Position Velocity report and may be useful for future attitude or in-track/cross-track analysis.
Sample run inputs used for this dataset:
- Start day number: `170.67`
- End day number: `170.78`
- Time step: `0.00073` days
- Face: `0`
- Tilt angle: `0` degrees
