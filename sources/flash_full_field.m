function [ ] = flash_full_field( input_args )
%FLASH_FULL_FIELD Summary of this function goes here
%   Detailed explanation goes here
    
    [~, mname, ~] = fileparts(mfilename('fullpath'));
    
    % 1. Parse input parameters.
    % % Define default value for each input parameter.
    flash_duration = 0.5; % sec
    post_flash_duration = 0.5; % sec
    initial_adaptation_duration = 0.0; % sec
    trial_adaptation_duration = 1.5; % sec
    background_intensity = 0.1;
    nb_repetitions = 5 * 30;
    dmd_width = 1024; % px
    dmd_height = 768; % px
    dmd_frame_rate = 60.0; % Hz
    dmd_inversed_polarity = false;
    output_foldername = pwd;
    % % Define input parser.
    parser = inputParser;
    parser.addParameter('flash_duration', flash_duration);
    parser.addParameter('post_flash_duration', post_flash_duration);
    parser.addParameter('initial_adaptation_duration', initial_adaptation_duration);
    parser.addParameter('trial_adaptation_duration', trial_adaptation_duration);
    parser.addParameter('background_intensity', background_intensity);
    parser.addParameter('nb_repetitions', nb_repetitions);
    parser.addParameter('dmd_width', dmd_width);
    parser.addParameter('dmd_height', dmd_height);
    parser.addParameter('dmd_frame_rate', dmd_frame_rate);
    parser.addParameter('dmd_inversed_polarity', dmd_inversed_polarity);
    parser.addParameter('output_foldername', output_foldername);
    % % Parse input arguments.
    parser.parse(input_args{:});
    % % Retrieve values of input parameters.
    args = parser.Results;
    
    % ...
    n_flash = ceil(args.flash_duration * args.dmd_frame_rate);
    lum_flash = ones(n_flash, 1);
    n_post_flash = ceil(args.post_flash_duration * args.dmd_frame_rate);
    lum_post_flash = zeros(n_post_flash, 1);
    lum_stim = [
        lum_flash;
        lum_post_flash;
    ];
    % % Rescale stimulus luminance profile.
    if args.background_intensity < 1.0
        args.background_intensity = 256.0 * args.background_intensity;
    end
    disp(['Background intensity: ', num2str(args.background_intensity)]);
    lum_stim = lum_stim * (256.0 - args.background_intensity);
    lum_stim = lum_stim + args.background_intensity;
    lum_stim = uint8(lum_stim - 0.5);
    
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
    lum_trial = args.background_intensity * ones(n_trial, 1);
    lum_trial = uint8(lum_trial - 0.5);
    
    % Generate full-field profile.
    profiles = cell(1 + 2 * nb_repetitions, 1);
    profiles{1} = lum_init;
    for repetition_id = 1:nb_repetitions
        profiles{2 * repetition_id + 0} = lum_trial;
        profiles{2 * repetition_id + 1} = lum_stim;
    end
    lum = cell2mat(profiles);
    
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
    repetitions = zeros(nb_repetitions, 3);
    nbs_frames = cellfun(@(lum) length(lum), profiles);
    frame_ids = cumsum(nbs_frames);
    for repetition_id = 1:nb_repetitions
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

