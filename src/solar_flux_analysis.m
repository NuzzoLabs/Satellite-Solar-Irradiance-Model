function results = solar_flux_analysis(inputFile, n_i, n_f, step, face, sigma, makePlots)
%SOLAR_FLUX_ANALYSIS Estimate incident solar flux on a satellite surface.
%
% results = solar_flux_analysis(inputFile, n_i, n_f, step, face, sigma, makePlots)
%
% Inputs
%   inputFile  - path to STK Fixed Position Velocity report exported as TSV/CSV
%   n_i        - start day number of year, including UTC fraction if needed
%   n_f        - end day number of year, including UTC fraction if needed
%   step       - time step in days
%   face       - face identifier retained for future body-frame attitude logic
%   sigma      - tilt angle relative to local Earth surface, degrees
%   makePlots  - true/false flag for plotting
%
% Assumed input columns
%   1 date, 2 time, 3 ECEF X km, 4 ECEF Y km, 5 ECEF Z km,
%   6 ECEF VX km/s, 7 ECEF VY km/s, 8 ECEF VZ km/s
%
% Notes
%   This code was refactored from a 2019 MATLAB/STK project. It assumes a
%   spherical Earth and a nadir-pointing attitude approximation.

arguments
    inputFile (1,:) char
    n_i (1,1) double
    n_f (1,1) double
    step (1,1) double
    face (1,1) double = 0 %#ok<INUSA>
    sigma (1,1) double = 0
    makePlots (1,1) logical = true
end

RE = 6371;              % Earth radius, km
SF = 1377;              % Solar flux at 1 AU, W/m^2
e_earth = 0.0174;       % Earth orbital eccentricity
theta_rate = 0.9863;    % Earth angular motion around sun, deg/day
rpe = 147.1e6;          % Earth perihelion distance, km
AU_km = 149597871;      % Astronomical unit, km

% Build day-number vector.
n = build_day_vector(n_i, n_f, step);

% Read STK data. MATLAB readmatrix ignores nonnumeric date/time columns when
% NumHeaderLines is not needed for a clean mixed TSV, so table import is safer.
opts = detectImportOptions(inputFile, 'FileType', 'text', 'Delimiter', '\t');
opts.VariableNames = {'Date','Time','ECEF_X_km','ECEF_Y_km','ECEF_Z_km','ECEF_VX_km_s','ECEF_VY_km_s','ECEF_VZ_km_s'};
T = readtable(inputFile, opts);

R_ECEF_x = T.ECEF_X_km.';
R_ECEF_y = T.ECEF_Y_km.';
R_ECEF_z = T.ECEF_Z_km.';

% Align time vector and STK rows if small rounding mismatch exists.
N = min([numel(n), numel(R_ECEF_x)]);
n = n(1:N);
R_ECEF_x = R_ECEF_x(1:N);
R_ECEF_y = R_ECEF_y(1:N);
R_ECEF_z = R_ECEF_z(1:N);

RMAG = sqrt(R_ECEF_x.^2 + R_ECEF_y.^2 + R_ECEF_z.^2);
ALT = RMAG - RE;

% Longitude using atan2d handles all quadrants.
longitude_deg = atan2d(R_ECEF_y, R_ECEF_x);

% Latitude assuming spherical Earth.
recef_unit_x = R_ECEF_x ./ RMAG;
recef_unit_y = R_ECEF_y ./ RMAG;
recef_unit_z = R_ECEF_z ./ RMAG;
Rsurf_x = recef_unit_x * RE;
Rsurf_y = recef_unit_y * RE;
Rsurf_z = recef_unit_z * RE;
rsurf_plane_mag = sqrt(Rsurf_x.^2 + Rsurf_y.^2);
latitude_deg = atan2d(Rsurf_z, rsurf_plane_mag);

% Solar declination.
Soldec = 23.45 * sind((360/365) * (n - 81));

% Earth true anomaly and space solar flux at Earth orbit.
theta = build_true_anomaly_vector(n_i, n_f, step, theta_rate);
theta = theta(1:min(numel(theta), N));
n = n(1:numel(theta));
longitude_deg = longitude_deg(1:numel(theta));
latitude_deg = latitude_deg(1:numel(theta));
RMAG = RMAG(1:numel(theta));
ALT = ALT(1:numel(theta));
Soldec = Soldec(1:numel(theta));

C = rpe * (1 + e_earth * cosd(0));
r_E = C ./ (1 + e_earth * cosd(theta));
Solflux_sp = SF * (AU_km ./ r_E).^2;

% Time-of-day vector and hour angles.
time_hr = build_time_hours(n_i, n_f, step, n);
Hourangle_pm = (12 - time_hr) * 15;
HRsat = satellite_hour_angle(Hourangle_pm, longitude_deg);

% Solar and geometric angles.
Beta = asind(cosd(latitude_deg) .* cosd(Soldec) .* cosd(HRsat) + sind(latitude_deg) .* sind(Soldec));
Dh = sqrt(RMAG.^2 - RE^2);
gamma = asind(Dh ./ RMAG);
l = sqrt(2*RE^2 - 2*RE^2 .* cosd(gamma));
s = acosd((l.^2 - ALT.^2 - Dh.^2) ./ (-2 .* ALT .* Dh));
Beta_corrected = Beta + (90 - s);
Phis = asind((cosd(Soldec) .* sind(HRsat)) ./ cosd(Beta));

A = diff(longitude_deg);
B = diff(latitude_deg);
phic = atan2d(B, A);
Phic = [phic, phic(end)];
DIF = Phis - Phic;

theta_inc = acosd(cosd(Beta) .* cosd(Phis - Phic) .* sind(sigma) + sind(Beta) .* cosd(sigma));

SOLF = zeros(size(n));
if sigma == 0
    for q = 1:numel(n)
        if HRsat(q) >= 90 || HRsat(q) <= -90
            SOLF(q) = 0;
        else
            SOLF(q) = Solflux_sp(q) * cosd(theta_inc(q));
        end
    end
elseif sigma == 90
    for q = 1:numel(n)
        if Beta_corrected(q) <= 0
            SOLF(q) = 0;
        else
            SOLF(q) = Solflux_sp(q) * cosd(theta_inc(q));
        end
    end
else
    for q = 1:numel(n)
        if Beta_corrected(q) <= 0
            SOLF(q) = 0;
        else
            SOLF(q) = max(0, Solflux_sp(q) * cosd(theta_inc(q)));
        end
    end
end

results = struct();
results.dayNumber = n;
results.longitude_deg = longitude_deg;
results.latitude_deg = latitude_deg;
results.altitude_km = ALT;
results.solarFluxSpace_W_m2 = Solflux_sp;
results.solarFluxSurface_W_m2 = SOLF;
results.solarDeclination_deg = Soldec;
results.satelliteHourAngle_deg = HRsat;
results.solarAltitude_deg = Beta;
results.solarAzimuth_deg = Phis;
results.groundTrackSlope_deg = Phic;
results.incidenceAngle_deg = theta_inc;
results.azimuthDifference_deg = DIF;

if makePlots
    plot_results(results);
end
end

function n = build_day_vector(n_i, n_f, step)
if n_f > n_i
    n = n_i:step:n_f;
else
    mat = n_i:step:366;
    if mat(end) == 366
        n1 = n_i:step:(366-step);
        n2 = 1:step:n_f;
    else
        n1 = mat;
        n2 = (step - (366 - mat(end))) + 1:step:n_f;
    end
    n = [n1, n2];
end
end

function theta = build_true_anomaly_vector(n_i, n_f, step, theta_rate)
if n_i >= 3
    theta_i = (n_i - 3) * theta_rate;
else
    theta_i = 360 - (3 - n_i) * theta_rate;
end

if n_f >= 3
    theta_f = (n_f - 3) * theta_rate;
else
    theta_f = 360 - (3 - n_f) * theta_rate;
end

theta_step = theta_rate * step;
if n_i >= 3 && theta_f >= theta_i
    theta = theta_i:theta_step:theta_f;
else
    t1 = theta_i:theta_step:(360 - theta_step);
    t2 = (theta_step - (360 - t1(end))):theta_step:theta_f;
    theta = [t1, t2];
end
end

function time_hr = build_time_hours(n_i, n_f, step, n)
t_start = (n_i - floor(n_i)) * 24;
t_fin = (n_f - floor(n_f)) * 24;
step_hr = 24 * step;

if floor(n_f) == 1 + floor(n_i)
    mat2 = t_start:step_hr:24;
    if mat2(end) == 24
        timei = t_start:step_hr:(24-step_hr);
        timef = 0:step_hr:t_fin;
    else
        timei = mat2;
        timef = step_hr - (24 - timei(end)):step_hr:t_fin;
    end
    time_hr = [timei, timef];
elseif floor(n_f) == floor(n_i)
    time_hr = t_start:step_hr:t_fin;
else
    time_hr = zeros(1, numel(n));
    t = floor(n);
    for m = 1:numel(t)
        if m == 1
            time_hr(m) = t_start;
        elseif m == numel(t)
            time_hr(m) = t_fin;
        elseif (t(m) == t(m-1))
            time_hr(m) = time_hr(m-1) + step_hr;
        else
            time_hr(m) = step_hr - (24 - time_hr(m-1));
        end
    end
end

time_hr = time_hr(1:min(numel(time_hr), numel(n)));
end

function HRsat = satellite_hour_angle(Hourangle_pm, longitude_deg)
HRsat = zeros(size(Hourangle_pm));
for p = 1:numel(Hourangle_pm)
    lam = longitude_deg(p);
    hapm = Hourangle_pm(p);
    if hapm < 0 && lam > 0 && ((lam + abs(hapm)) <= 180)
        HRsat(p) = -lam + hapm;
    elseif hapm < 0 && lam > 0 && ((lam + abs(hapm)) > 180)
        HRsat(p) = 360 - (abs(hapm) + lam);
    elseif hapm < 0 && lam < 0
        HRsat(p) = -lam + hapm;
    elseif hapm > 0 && lam > 0
        HRsat(p) = -lam + hapm;
    elseif hapm > 0 && lam < 0 && ((abs(lam) + hapm) <= 180)
        HRsat(p) = -lam + hapm;
    elseif hapm > 0 && lam < 0 && ((abs(lam) + hapm) > 180)
        HRsat(p) = 360 - (abs(lam) + hapm);
    else
        HRsat(p) = -lam + hapm;
    end
end
end

function plot_results(results)
figure; plot(results.longitude_deg, results.latitude_deg, 'LineWidth', 2); grid on;
title('Ground Track of Orbit'); xlabel('Longitude (deg)'); ylabel('Latitude (deg)');

figure; plot(results.dayNumber, results.solarFluxSpace_W_m2); grid on;
title('Solar Flux in Space at Earth Orbit'); xlabel('Day Number'); ylabel('Solar Flux (W/m^2)');

figure; plot(results.dayNumber, results.solarFluxSurface_W_m2); grid on;
title('Solar Flux on Satellite Surface as Function of Time'); xlabel('Day Number'); ylabel('Solar Flux (W/m^2)');

figure; plot(results.longitude_deg, results.solarFluxSurface_W_m2); grid on;
title('Solar Flux on Satellite Surface as Function of Longitude'); xlabel('Longitude (deg)'); ylabel('Solar Flux (W/m^2)');

figure; plot(results.latitude_deg, results.solarFluxSurface_W_m2); grid on;
title('Solar Flux on Satellite Surface as Function of Latitude'); xlabel('Latitude (deg)'); ylabel('Solar Flux (W/m^2)');

figure; plot(results.dayNumber, results.longitude_deg, results.dayNumber, results.latitude_deg); grid on;
title('Latitude and Longitude as a Function of Time'); xlabel('Day Number'); ylabel('Lat and Long (deg)'); legend('Longitude','Latitude');

figure; plot(results.dayNumber, results.satelliteHourAngle_deg); grid on;
title('Satellite Hour Angle as Function of Time'); xlabel('Day Number'); ylabel('Hour Angle (deg)');

figure; plot(results.dayNumber, results.incidenceAngle_deg); grid on;
title('Solar Incidence Angle on Surface of Satellite'); xlabel('Day Number'); ylabel('Incidence Angle (deg)');
end
