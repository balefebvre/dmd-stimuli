function [ ] = binary_time_series( input_args )
% Function for generating a DMD stimulus based on
% a binary input X (loaded from a .mat file).

[~, mname, ~] = fileparts(mfilename('fullpath'));

% 1. Parse input parameters.
% % Define default value for each input parameter.
on_duration = 0.04; % sec
post_on_duration = 0.08; % sec
off_duration = 0.04; % sec
post_off_duration = 0.08; % sec
initial_adaptation_duration = 0.0; % sec
trial_adaptation_duration = 0.0; % sec
background_intensity = 0.1;
total_duration = 3600; % sec
nb_repetitions = 1;
dmd_width = 1024; % px
dmd_height = 768; % px
dmd_frame_rate = 100.0; % Hz
dmd_inversed_polarity = true;
input_path = '';
output_foldername = pwd;
% % Define input parser.
parser = inputParser;
parser.addParameter('on_duration', on_duration);
parser.addParameter('post_on_duration', post_on_duration);
parser.addParameter('off_duration', off_duration);
parser.addParameter('post_off_duration', post_off_duration);
parser.addParameter('initial_adaptation_duration', initial_adaptation_duration);
parser.addParameter('trial_adaptation_duration', trial_adaptation_duration);
parser.addParameter('background_intensity', background_intensity);
parser.addParameter('total_duration', total_duration);
parser.addParameter('nb_repetitions', nb_repetitions);
parser.addParameter('dmd_width', dmd_width);
parser.addParameter('dmd_height', dmd_height);
parser.addParameter('dmd_frame_rate', dmd_frame_rate);
parser.addParameter('dmd_inversed_polarity', dmd_inversed_polarity);
parser.addParameter('input_path', input_path);
parser.addParameter('output_foldername', output_foldername);
% % Parse input arguments.
parser.parse(input_args{:});
% % Retrieve values of input parameters.
args = parser.Results;

% Defining luminance levels.
% on signal (1) in input
n_on = ceil(args.on_duration * args.dmd_frame_rate);
lum_on = ones(n_on, 1);
n_post_on = ceil(args.post_on_duration * args.dmd_frame_rate);
lum_post_on = zeros(n_post_on, 1);
lum_stim = [
    lum_on;
    lum_post_on;
];
% off signal (0) in input
n_off = ceil(args.off_duration * args.dmd_frame_rate);
lum_off = zeros(n_off, 1);
n_post_off = ceil(args.post_off_duration * args.dmd_frame_rate);
lum_post_off = ones(n_post_off, 1);
lum_sil = [
    lum_off;
    lum_post_off;
];

% % Rescale stimulus luminance profile.
if args.background_intensity < 1.0
    args.background_intensity = 256.0 * args.background_intensity;
end
disp(['Background intensity: ', num2str(args.background_intensity)]);

lum_stim = lum_stim * (256.0 - args.background_intensity);
lum_stim = lum_stim + args.background_intensity;
lum_stim = uint8(lum_stim - 0.5);

lum_sil = lum_sil * (256.0 - args.background_intensity);
lum_sil = lum_sil + args.background_intensity;
lum_sil = uint8(lum_sil - 0.5);

% Compute luminance statistics.
disp(['Minimal luminance: ', num2str(min(lum_stim))]);
disp(['Median luminance: ', num2str(median(lum_stim))]);
disp(['Maximal luminance: ', num2str(max(lum_stim))]);

% Compute stimulus duration.
n_stim = length(lum_stim);
d_stim = n_stim / args.dmd_frame_rate;
disp('Stimulus duration:');
disp(['  ', num2str(d_stim), ' sec']);
if d_stim > 60.0
    dm = d_stim / 60.0; % minutes
    ds = d_stim - dm * 60.0; %s econds
    disp(['  ', dm, ' min ', ds, ' sec']);
end

% Make output directory.
if ~exist(output_foldername, 'dir')
    mkdir(output_foldername);
end

% Generate full-field profile used at the beginning.
n_init = ceil(args.initial_adaptation_duration * args.dmd_frame_rate);
lum_init = args.background_intensity * ones(n_init, 1);
lum_init = uint8(lum_init - 0.5);

% Generate full-field profile used before each repetition.
n_trial = ceil(args.trial_adaptation_duration * args.dmd_frame_rate);
if n_trial
    lum_trial = args.background_intensity * ones(n_trial, 1);
    lum_trial = uint8(lum_trial - 0.5);
end

% Generate full-field profile.
nb_rep = args.nb_repetitions;
profiles = cell(1 + 2 * nb_rep, 1);
profiles{1} = lum_init;

N = args.total_duration / d_stim; % given that 0 and 1 last the same (on+post_on)

% Load time-series X, containing 0s and 1s.
try
    load(args.input_path, 'X')
    disp('*** Input file successfully loaded. ***')
catch
    error('*** Input file could not be loaded. Check if the path is correct (file extension included) and if the file exists. ***')
end

% Adjust the input stimulus if necessary.
if(numel(X) > N)
    X = X(1:N); 
    warning('Input X had to be clipped to fir N steps.')
end

% Replace zeros with lum_sil and ones with lum_stim.
for repetition_id = 1:nb_repetitions 
    for i = 1:numel(X)
        if(X(i))
            profiles{i} = lum_stim;
        else
            profiles{i} = lum_sil;
        end
    end
end
lum = cell2mat(profiles);

% Display an example of the stimulus that will be saved.
show_trace(lum, 1000, args.dmd_frame_rate);

% Compute and display stimulus duration.
nb_frames = length(lum);
d_total =  nb_frames / args.dmd_frame_rate;
disp('Total duration:');
disp(['  ', num2str(d_total), ' sec']);
if d_total > 60.0
    dm = floor(d_total / 60.0); % minutes
    ds = d_total - dm * 60.0; % seconds
    disp(['  ', num2str(dm), ' min ', num2str(ds), ' sec']);
end

% Write .bin file.
% % Open .bin file.
bin_filename = [mname, '_', num2str(args.dmd_frame_rate), 'hz.bin'];
bin_pathname = fullfile(args.output_foldername, bin_filename);
permission = 'w'; % writing mode, discard existing contents
machine_format = 'l'; % IEEE floating point with little-endian byte ordering
bin_fid = fopen(bin_pathname, permission, machine_format);
% % Write .bin file header.
nb_images = 256; % number of images
nb_bits = 8; % number of bits
bin_header = [args.dmd_width, args.dmd_height, nb_images, nb_bits];
fwrite(bin_fid, bin_header, 'int16');
% % Write .bin file images.
bin_image = zeros(args.dmd_height, args.dmd_width, 'uint8');
for i = 1:nb_images
    if args.dmd_inversed_polarity
        fwrite(bin_fid, 255 - bin_image(:), 'uint8');
    else
        fwrite(bin_fid, bin_image(:), 'uint8');
    end
    bin_image = bin_image + 1;
end
% % Close .bin file.
fclose(bin_fid);

% Write .vec file.
% % Open .vec file.
vec_filename = [mname, '_', num2str(args.dmd_frame_rate), 'hz.vec'];
vec_pathname = fullfile(args.output_foldername, vec_filename);
permission = 'w'; % writing mode, discard existing contents
machine_format = 'l'; % IEEE floating point with little-endian byte ordering
vec_fid = fopen(vec_pathname, permission, machine_format);
% % Write .vec file header.
vec_header = [0, nb_frames, 0, 0, 0];
fprintf(vec_fid, '%g %g %g %g %g\n', vec_header);
% % Write .vec file image indices.
for frame_id = 1:nb_frames
    image_id = lum(frame_id);
    vec_frame = [0, image_id, 0, 0, 0];
    fprintf(vec_fid, '%g %g %g %g %g\n', vec_frame);
end
% % Close .vec file.
fclose(vec_fid);

% Define repetitions data.
repetitions = zeros(nb_rep, 3);
nbs_frames = cellfun(@(lum) length(lum), profiles);
frame_ids = cumsum(nbs_frames);
for repetition_id = 1:nb_rep
    repetitions(repetition_id, 1) = repetition_id;
    start_frame_id = frame_ids(2 * repetition_id) + 1;
    repetitions(repetition_id, 2) = start_frame_id;
    end_frame_id = frame_ids(2 * repetition_id + 1);
    repetitions(repetition_id, 3) = end_frame_id;
end

% Write repetitions file.
% % Open repetitions file.
rep_filename = [mname, '_repetitions.csv'];
rep_pathname = fullfile(args.output_foldername, rep_filename);
rep_fid = fopen(rep_pathname, 'w');
% % Write repetitions file header.
rep_header = 'repetitionId;startFrameId;endFrameId';
fprintf(rep_fid, '%s\r\n', rep_header);
% % Close repetitions file.
fclose(rep_fid);
% % Write repetitions file data.
dlmwrite(rep_pathname, repetitions, '-append', 'delimiter', ';', 'newline', 'pc');

% Define stimulus data.
nb_stim_frames = length(lum_stim);
stimulus = zeros(nb_stim_frames, 2);
stimulus(:, 1) = 1:nb_stim_frames;
stimulus(:, 2) = lum_stim;

% Write stimulus file.
% % Open stimulus file.
stim_filename = [mname, '_stimulus.csv'];
stim_pathname = fullfile(args.output_foldername, stim_filename);
stim_fid = fopen(stim_pathname, 'w');
% % Write stimulus file header.
stim_header = 'frameId;luminance';
fprintf(stim_fid, '%s\r\n', stim_header);
% % Close stimulus file.
fclose(stim_fid);
% % Write stimulus file data.
dlmwrite(stim_pathname, stimulus, '-append', 'delimiter', ';', 'newline', 'pc');

% Define parameter data.
args.stimulus_duration = d_stim;
args.total_duration = d_total;
args.bin_filename = bin_filename;
args.vec_filename = vec_filename;
args.repetition_filename = rep_filename;
args.stimulus_filename = stim_filename;
% % Discard unused warning.
args; %#ok

% Write parameters file.
mat_filename = fullfile(args.output_foldername, [mname, '_parameters.mat']);
save(mat_filename, '-struct', 'args');

return

end

