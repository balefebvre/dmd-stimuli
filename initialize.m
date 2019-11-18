% Initialize the output folder and add paths.

% path to where /dmd-stimuli folder is located
[mfoldername, ~, ~] = fileparts(mfilename('fullpath'));

% creating an output folder in the same path
outputs_foldername = fullfile(mfoldername, 'outputs');
if ~exist(outputs_foldername, 'dir')
    mkdir(outputs_foldername);
end

addpath(genpath(mfoldername))
% addpath(fullfile(mfoldername, 'sources'));
% addpath(fullfile(mfoldername, 'sources', 'euler_full_field'));
% addpath(fullfile(mfoldername, 'sources', 'moving_bars'));
% addpath(fullfile(mfoldername, 'examples'));