[mfoldername, ~, ~] = fileparts(mfilename('fullpath'));

input_args = {...
    'flash_duration', 0.5, ... % sec
    'post_flash_duration', 0.5, ... % sec
    'background_intensity', 0.1, ...
    'initial_adaptation_duration', 0.0, ... % sec
    'trial_adaptation_duration', 1.0, ... % sec
    'nb_repetitions', 5 * 30, ...
    'dmd_width', 1024, ... % px
    'dmd_height', 768, ... % px
    'dmd_frame_rate', 100.0, ... % Hz
    'dmd_inversed_polarity', false, ...
    'output_foldername', fullfile(pwd, 'outputs'), ...
};

flash_full_field(input_args);