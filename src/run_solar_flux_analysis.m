%% Run Satellite Solar Flux Analysis
% Interactive runner for the included STK sample dataset.

clc;
clear;

inputFile = fullfile('data', 'orbit_data_no_maneuver.tsv');

y1 = ['What is the day number at start of STK simulation? ', ...
      '(Use UTC time as fraction of day): '];
n_i = input(y1);

y2 = 'What is the day number at the end of the STK simulation data set? ';
n_f = input(y2);

y3 = ['What is the time step for the STK data in number of days? ', ...
      'Example: 0.00073: '];
step = input(y3);

y4 = 'What is the face that you want to analyze? ';
face = input(y4); %#ok<NASGU> retained for future face-specific attitude logic

y5 = ['What is the tilt angle relative to the surface of Earth? ', ...
      'Input initial tilt angle if attitude is time-dependent: '];
sigma = input(y5);

results = solar_flux_analysis(inputFile, n_i, n_f, step, face, sigma, true);

fprintf('\nAnalysis complete. Computed %d solar-flux samples.\n', numel(results.solarFluxSurface_W_m2));
