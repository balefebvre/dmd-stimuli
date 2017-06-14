[mfoldername, ~, ~] = fileparts(mfilename('fullpath'));

outputs_foldername = fullfile(mfoldername, 'outputs');
if ~exist(outputs_foldername, 'dir')
    mkdir(outputs_foldername);
end

addpath(fullfile(mfoldername, 'sources'));
addpath(fullfile(mfoldername, 'sources', 'euler_full_field'));
addpath(fullfile(mfoldername, 'sources', 'moving_bars'));
addpath(fullfile(mfoldername, 'examples'));