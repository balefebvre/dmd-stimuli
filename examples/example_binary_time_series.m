% Path setup.
[mfoldername, ~, ~] = fileparts(mfilename('fullpath'));

% Set input argument values.
input_args = {...
    'on_duration', 0.04, ... % sec % 40 ms of ones
    'off_duration', 0.12, ... % sec % 120 ms of zeros
    'post_on_duration', 0.08, ... % sec % 80 ms of zeros after ones
    'post_off_duration', 0.0, ... % sec % 0 ms of ones after zeros
    'background_intensity', 0, ... %
    'initial_adaptation_duration', 0.0, ... % sec
    'trial_adaptation_duration', 0.0, ... % sec
    'total_duration', 3600, ... % sec
    'nb_repetitions', 1, ...
    'dmd_width', 1024, ... % px
    'dmd_height', 768, ... % px
    'dmd_frame_rate', 100.0, ... % Hz
    'dmd_inversed_polarity', true, ...
    'input_path', 'C:\Users\dmd-stimuli\input.mat'...
    'output_foldername', fullfile(pwd, 'outputs'), ...
};

% Generate stimulus.
binary_time_series(input_args);

% Show the generated stimulus. 
dmd_frame_rate = 100; 
fprintf('Displaying generated stimulus at %d Hz.\n', dmd_frame_rate)
path = strcat(input_args{end}, '\binary_time_series_', num2str(dmd_frame_rate), 'hz');
play_stim(path, dmd_frame_rate)
