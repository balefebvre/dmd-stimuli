[mfoldername, ~, ~] = fileparts(mfilename('fullpath'));

addpath(fullfile(mfoldername, 'sources'));
addpath(fullfile(mfoldername, 'sources', 'euler_full_field'));
addpath(fullfile(mfoldername, 'examples'));