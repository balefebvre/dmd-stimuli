% Path setup.
[mfoldername, ~, ~] = fileparts(mfilename('fullpath'));

% Set input argument values.
input_args = {...
    'rectangle_width', 600, ...
    'rectangle_height', 270, ...
    'rectangle_duration', 0.3, ... % sec
    'post_rect_duration', 0.3, ... % sec
    'initial_adaptation_duration', 0.0, ... % sec 
    'trial_adaptation_duration', 0.0, ... % sec
    'rectangle_intensity', 0.9, ... % pixel intensity (1 to 255 or 0.0 to 1.0)
    'background_intensity', 0.0, ... % 
    'nb_repetitions', 3, ...
    'dmd_width', 1024, ... % px
    'dmd_height', 768, ... % px
    'dmd_frame_rate', 100.0, ... % Hz
    'dmd_inversed_polarity', false, ...
    'output_foldername', fullfile(pwd, 'outputs'), ...
};

% Generate stimulus.
rectangle_central_field(input_args);

% Show the generated stimulus. 
dmd_frame_rate = 100; 
fprintf('Displaying generated stimulus at %d Hz.\n', dmd_frame_rate)
path = strcat(input_args{end}, '\rectangle_central_field_', ...
                                num2str(dmd_frame_rate), 'hz');
play_stim(path, dmd_frame_rate)

