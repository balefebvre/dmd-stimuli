function [ ] = natural_full_field( input_args )
%NATURAL_FULL_FIELD Generate a natural full-field DMD stimulus.
%   Detailed explanation goes here
    
    [mfoldername, mname, ~] = fileparts(mfilename('fullpath'));
    
    % 1. Parse input parameters.

    % % Define default value for each input parameter.
    input_foldername = fullfile(mfoldername, '..', 'data');
    input_filename = 'green_nature_waterfall.jpg';
    background_intensity = 0.3;
    initial_adaptation_duration = 20.0; % sec
    trial_adaptation_duration = 10.0; % sec
    nb_repetitions = 10;
    dmd_width = 1024; % px
    dmd_height = 768; % px
    dmd_frame_rate = 60; % Hz
    output_foldername = pwd;
    
    % % Define input parser.
    parser = inputParser;
    parser.addParameter('input_foldername', input_foldername);
    parser.addParameter('input_filename', input_filename);
    parser.addParameter('background_intensity', background_intensity);
    parser.addParameter('initial_adaptation_duration', initial_adaptation_duration);
    parser.addParameter('trial_adaptation_duration', trial_adaptation_duration);
    parser.addParameter('nb_repetitions', nb_repetitions);
    parser.addParameter('dmd_width', dmd_width);
    parser.addParameter('dmd_height', dmd_height);
    parser.addParameter('dmd_frame_rate', dmd_frame_rate);
    parser.addParameter('output_foldername', output_foldername);
    
    % % Parse input arguments.
    parser.parse(input_args{:});
    
    % % Retrieve value for each input parameter.
    results = parser.Results;
    input_foldername = results.input_foldername;
    input_filename = results.input_filename;
    background_intensity = results.background_intensity;
    initial_adaptation_duration = results.initial_adaptation_duration;
    trial_adaptation_duration = results.trial_adaptation_duration;
    nb_repetitions = results.nb_repetitions;
    dmd_width = results.dmd_width;
    dmd_height = results.dmd_height;
    dmd_frame_rate = results.dmd_frame_rate;
    output_foldername = results.output_foldername;
    
    
    
    % Define path of graphics file.
    path = fullfile(input_foldername, input_filename);
    disp(['Image path: ', path]);
    
    % Read image from graphics file.
    A = imread(path);
    
    % Get size of image.
    h = size(A, 1); % image height
    disp(['Image height: ', num2str(h)]);
    w = size(A, 2); % image width
    disp(['Image width: ', num2str(w)]);
    s = [h, w]; % image size
    
    % Define interpolation points.
    n = 5000;
    t = linspace(0.0, 1.0, n);
    p_start = [0.11, 0.1] .* s;
    p_end = [0.9, 0.9] .* s;
    p = repmat(p_start, n, 1) + repmat(p_end - p_start, n, 1) .* repmat(t', 1, 2);
    
    % Make output directory.
    mkdir(output_foldername);
    
    % Plot first control figure.
    figure('visible', 'off');
    hold('on');
    image(A);
    plot(p(:, 2), p(:, 1), 'w-');
    xlim([0, w]);
    ylim([0, h]);
    set(gca, 'Ydir', 'reverse');
    hold('off');
    xlabel('pixel');
    ylabel('pixel');
    saveas(gcf, fullfile(output_foldername, [mname, '_control_1.svg']));
    close();
    
    % Get interpolated luminance values.
    [I, J] = meshgrid(1:w, 1:h);
    L = mean(A, 3);
    i = p(:, 2);
    j = p(:, 1);
    l_base = interp2(I, J, L, i, j);
    
    % Normalize luminance values.
    if background_intensity < 1.0
        background_intensity = 256.0 * background_intensity;
    else
        background_intensity = 1.0 * background_intensity;
    end
    l_base = l_base + (background_intensity - min(l_base));
    l_base = l_base * (256.0 / max(l_base));
    l_base = uint8(l_base);
    disp(['Maximum light intensity: ', num2str(max(l_base))]);
    disp(['Median light intensity: ', num2str(median(l_base))]);
    disp(['Minimum light intensity: ', num2str(min(l_base))]);
    
    % Plot second control figure.
    x = (1:n) / dmd_frame_rate;
    y = l_base;
    figure('visible', 'off');
    hold('on');
    plot(x, y, 'k-');
    hold('off');
    xlim([min(x) max(x)]);
    ylim([min(y) max(y)]);
    xlabel('time (sec)');
    ylabel('luminance (arb. unit)');
    saveas(gcf, fullfile(output_foldername, [mname, '_control_2.svg']));
    close();
    
    % Define adaptation durations.
    d_init = initial_adaptation_duration;
    d_trial = trial_adaptation_duration;
    n_init = floor(d_init * dmd_frame_rate);
    n_trial = floor(d_trial * dmd_frame_rate);
    d_init = n_init / dmd_frame_rate;
    d_trial = n_trial / dmd_frame_rate;
    disp(['Adaptation duration (initial): ', num2str(d_init), ' sec']);
    disp(['Adaptation duration (before repetition): ', num2str(d_trial), ' sec']);
    
    % Define luminance profiles for retina adaptation.
    median_intensity = median(l_base);
    l_init = median_intensity * ones(n_init, 1, 'uint8');
    l_trial = median_intensity * ones(n_trial, 1, 'uint8');
    
    % Define apparition order of luminance profiles.
    nb_traces = 1 + 2 * nb_repetitions;
    traces = cell(nb_traces, 1);
    traces{1} = l_init;
    for i = 1:nb_repetitions
        traces{2 * i + 0} = l_trial;
        traces{2 * i + 1} = l_base;
    end
    
    % Compute and display stimulus duration.
    l = cell2mat(traces);
    nb_frames = size(l, 1);
    d =  nb_frames / dmd_frame_rate;
    disp(['Duration: ', num2str(d), ' sec']);
    if d > 60.0
        dm = floor(d / 60.0);
        ds = d - dm * 60.0;
        disp(['          ', num2str(dm), ' min ', num2str(ds), ' sec']);
    end
    
    % Write .bin file.
    % % Open .bin file.
    bin_filename = [mname, '.bin'];
    permission = 'w'; % writing mode, discard existing contents
    machine_format = 'l'; % IEEE floating point with little-endian byte ordering
    bin_fid = fopen(fullfile(output_foldername, bin_filename), permission, machine_format);
    % % Write .bin file header.
    nb_images = 256; % number of images
    nb_bits = 8; % number of bits
    bin_header = [dmd_width, dmd_height, nb_images, nb_bits];
    fwrite(bin_fid, bin_header, 'int16');
    % % Write .bin file images.
    bin_image = zeros(dmd_height, dmd_width, 'uint8');
    for i = 1:nb_images
        fwrite(bin_fid, bin_image(:), 'uint8');
        bin_image = bin_image + 1;
    end
    % % Close .bin file.
    fclose(bin_fid);
    
    % Write .vec file.
    % % Open .vec file.
    vec_filename = [mname, '.vec'];
    permission = 'w'; % writing mode, discard existing contents
    machine_format = 'l'; % IEEE floating point with little-endian byte ordering
    vec_fid = fopen(fullfile(output_foldername, vec_filename), permission, machine_format);
    % % Write .vec file header.
    vec_header = [0, nb_frames, 0, 0, 0];
    fprintf(vec_fid, '%g %g %g %g %g\n', vec_header);
    % % Write .vec file image indices.
    for frame_id = 1:nb_frames
        image_id = l(frame_id);
        vec_frame = [0, image_id, 0, 0, 0];
        fprintf(vec_fid, '%g %g %g %g %g\n', vec_frame);    
    end
    % % Close .vec file.
    fclose(vec_fid);
    
    % Define repetitions data.
    repetitions = zeros(nb_repetitions, 3);
    nbs_frames = cellfun(@(l) size(l, 1), traces);
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
    rep_filename = fullfile(output_foldername, [mname, '_repetitions.csv']);
    rep_fid = fopen(rep_filename, 'w');
    % % Write repetitions file header.
    rep_header = 'repetitionId,startFrameId,endFrameId';
    fprintf(rep_fid, '%s\r\n', rep_header);
    % % Close repetitions file.
    fclose(rep_fid);
    % % Write repetitions file data.
    dlmwrite(rep_filename, repetitions, '-append', 'delimiter', ',');
    
    % Define stimulus data.
    nb_base_frames = size(l_base, 1);
    stimulus = zeros(nb_base_frames, 2);
    stimulus(:, 1) = 1:nb_base_frames;
    stimulus(:, 2) = l_base;
    
    % Write stimulus file.
    % % Open stimulus file.
    stim_filename = fullfile(output_foldername, [mname, '_stimulus.csv']);
    stim_fid = fopen(stim_filename, 'w');
    % % Write stimulus file header.
    stim_header = 'frameId,luminance';
    fprintf(stim_fid, '%s\r\n', stim_header);
    % % Close stimulus file.
    fclose(stim_fid);
    % % Write stimulus file data.
    dlmwrite(stim_filename, stimulus, '-append', 'delimiter', ',');
    
end

