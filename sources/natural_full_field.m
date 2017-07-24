function [ ] = natural_full_field( input_args )
%NATURAL_FULL_FIELD Generate a natural full-field DMD stimulus.
%   Detailed explanation goes here
    
    [mfoldername, mname, ~] = fileparts(mfilename('fullpath'));
    
    % 1. Parse input parameters.
    % % Define input parser.
    parser = inputParser;
    parser.addParameter('input_foldername', fullfile(mfoldername, '..', 'data'));
    parser.addParameter('input_filename', 'green_nature_waterfall.jpg');
    parser.addParameter('stimulus_duration', 30.0); % sec
    parser.addParameter('background_intensity', 0.1);
    parser.addParameter('initial_adaptation_duration', 20.0); % sec
    parser.addParameter('trial_adaptation_duration', 10.0); % sec
    parser.addParameter('nb_repetitions', 10);
    parser.addParameter('dmd_width', 1024); % px
    parser.addParameter('dmd_height', 768); % px
    parser.addParameter('dmd_frame_rate', 60); % Hz
    parser.addParameter('dmd_inversed_polarity', false);
    parser.addParameter('output_foldername', pwd);
    % % Parse input arguments.
    parser.parse(input_args{:});
    % % Retrieve value of input parameters.
    args = parser.Results;
    
    
    % Define path of graphics file.
    path = fullfile(args.input_foldername, args.input_filename);
    disp(['Image path: ', path]);
    
    % Read image from graphics file.
    A = imread(path);
    
    % Get size of image.
    h = size(A, 1); % image height
    disp(['Image height: ', num2str(h)]);
    w = size(A, 2); % image width
    disp(['Image width: ', num2str(w)]);
    s = [h, w]; % image size
    
    % Define stimulus duration.
    d_stim = args.stimulus_duration;
    n_stim = ceil(d_stim * args.dmd_frame_rate);
    d_stim = n_stim / args.dmd_frame_rate;
    disp('Stimulus duration:');
    disp(['  ', num2str(d_stim), ' sec']);
    if d_stim > 60.0
        dm = floor(d_stim / 60.0); % minutes
        ds = d_stim - dm * 60.0; % seconds
        disp(['  ', num2str(dm), ' min ', num2str(ds), ' sec']);
    end
    
    % Define interpolation points.
    t = linspace(0.0, 1.0, n_stim);
    p_start = [0.11, 0.1] .* s;
    p_end = [0.9, 0.9] .* s;
    p = repmat(p_start, n_stim, 1) + repmat(p_end - p_start, n_stim, 1) .* repmat(t', 1, 2);
    
    % Make output directory.
    if ~exist(args.output_foldername, 'dir')
        mkdir(args.output_foldername);
    end
    
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
    saveas(gcf, fullfile(args.output_foldername, [mname, '_control_1.svg']));
    close();
    
    % Get interpolated luminance values.
    [I, J] = meshgrid(1:w, 1:h);
    L = mean(A, 3);
    i = p(:, 2);
    j = p(:, 1);
    l_base = interp2(I, J, L, i, j);
    
    % Normalize luminance values.
    if args.background_intensity < 1.0
        args.background_intensity = 256.0 * args.background_intensity;
    else
        args.background_intensity = 1.0 * args.background_intensity;
    end
    l_base = l_base - min(l_base);
    l_base = l_base / max(l_base);
    l_base = l_base * (256.0 - args.background_intensity);
    l_base = l_base + args.background_intensity;
    l_base = uint8(l_base - 0.5);
    
    maximum_intensity = max(l_base);
    median_intensity = median(l_base);
    minimum_intensity = min(l_base);
    disp(['Maximum light intensity: ', num2str(maximum_intensity)]);
    disp(['Median light intensity: ', num2str(median_intensity)]);
    disp(['Minimum light intensity: ', num2str(minimum_intensity)]);
    
    % Plot second control figure.
    x = (0:n_stim)  / args.dmd_frame_rate;
    y = [l_base; l_base(end)];
    figure('visible', 'off');
    hold('on');
    stairs(x, y, 'k-');
    hold('off');
    xlim([min(x), max(x)]);
    ylim([0.0-5.0, 255.0+5.0]);
    xlabel('time (sec)');
    ylabel('luminance (arb. unit)');
    saveas(gcf, fullfile(args.output_foldername, [mname, '_control_2.svg']));
    close();
    
    % Define adaptation durations.
    d_init = args.initial_adaptation_duration;
    d_trial = args.trial_adaptation_duration;
    n_init = floor(d_init * args.dmd_frame_rate);
    n_trial = floor(d_trial * args.dmd_frame_rate);
    d_init = n_init / args.dmd_frame_rate;
    d_trial = n_trial / args.dmd_frame_rate;
    disp(['Adaptation duration (initial): ', num2str(d_init), ' sec']);
    disp(['Adaptation duration (before repetition): ', num2str(d_trial), ' sec']);
    
    % Define luminance profiles for retina adaptation.
    lum_init = args.background_intensity * ones(n_init, 1);
    lum_init = uint8(lum_init - 0.5);
    l_trial = args.background_intensity * ones(n_trial, 1);
    l_trial = uint8(l_trial - 0.5);
    
    % Define apparition order of luminance profiles.
    nb_traces = 1 + 2 * args.nb_repetitions;
    traces = cell(nb_traces, 1);
    traces{1} = lum_init;
    for i = 1:args.nb_repetitions
        traces{2 * i + 0} = l_trial;
        traces{2 * i + 1} = l_base;
    end
    
    % Compute and display stimulus duration.
    l = cell2mat(traces);
    nb_frames = size(l, 1);
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
    bin_filename = [mname, '_', num2str(args.dmd_frame_rate),'hz.bin'];
    permission = 'w'; % writing mode, discard existing contents
    machine_format = 'l'; % IEEE floating point with little-endian byte ordering
    bin_fid = fopen(fullfile(args.output_foldername, bin_filename), permission, machine_format);
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
    permission = 'w'; % writing mode, discard existing contents
    machine_format = 'l'; % IEEE floating point with little-endian byte ordering
    vec_fid = fopen(fullfile(args.output_foldername, vec_filename), permission, machine_format);
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
    repetitions = zeros(args.nb_repetitions, 3);
    nbs_frames = cellfun(@(l) size(l, 1), traces);
    frame_ids = cumsum(nbs_frames);
    for repetition_id = 1:args.nb_repetitions
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
    rep_header = sprintf('repetitionId;startFrameId;endFrameId');
    fprintf(rep_fid, '%s\r\n', rep_header);
    % % Close repetitions file.
    fclose(rep_fid);
    % % Write repetitions file data.
    dlmwrite(rep_pathname, repetitions, '-append', 'delimiter', ';', 'newline', 'pc');
    
    % Define stimulus data.
    nb_base_frames = size(l_base, 1);
    stimulus = zeros(nb_base_frames, 2);
    stimulus(:, 1) = 1:nb_base_frames;
    stimulus(:, 2) = l_base;
    
    % Write stimulus file.
    % % Open stimulus file.
    stim_filename = [mname, '_stimulus.csv'];
    stim_pathname = fullfile(args.output_foldername, stim_filename);
    stim_fid = fopen(stim_pathname, 'w');
    % % Write stimulus file header.
    stim_header = sprintf('frameId;luminance');
    fprintf(stim_fid, '%s\r\n', stim_header);
    % % Close stimulus file.
    fclose(stim_fid);
    % % Write stimulus file data.
    dlmwrite(stim_pathname, stimulus, '-append', 'delimiter', ';', 'newline', 'pc');
    
    % Define parameter data.
    % % Internal parameters.
    args.minimum_intensity = minimum_intensity;
    args.median_intensity = median_intensity;
    args.maximum_intensity = maximum_intensity;
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
    
end

