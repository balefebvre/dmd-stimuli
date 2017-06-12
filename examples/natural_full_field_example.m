[mfoldername, ~, ~] = fileparts(mfilename('fullpath'));

input_args = {...
    'input_foldername', fullfile(mfoldername, '..', 'data'),...
    'input_filename', 'green_nature_waterfall.jpg',...
    'background_intensity', 0.3,...
    'initial_adaptation_duration', 20.0,... % sec
    'trial_adaptation_duration', 10.0,... % sec
    'nb_repetitions', 10,...
    'dmd_width', 1024,... % px
    'dmd_height', 768,... % px
    'dmd_frame_rate', 60.0, ... % Hz
    'output_foldername', fullfile(pwd, 'outputs'),...
};

natural_full_field(input_args);