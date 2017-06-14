function [ ] = moving_bars( input_args )
%MOVING_BARS Summary of this function goes here
%   Detailed explanation goes here
    
    [~, mname, ~] = fileparts(mfilename('fullpath'));
    
    % 1. Parse input parameters.
    % % Define input parser.
    parser = inputParser;
    parser.addParameter('seed', 42);
    parser.addParameter('nb_node_columns', 3);
    parser.addParameter('nb_node_rows', 3);
    parser.addParameter('node_spacing', 150.0); % µm
    parser.addParameter('trace_pad_length', 400.0); % µm
    parser.addParameter('bar_speed', 1000.0); % µm / sec
    parser.addParameter('bar_width', 300.0); % µm
    parser.addParameter('bar_length', 1000.0); % µm
    parser.addParameter('pixel_size', 2.3); % µm
    parser.addParameter('background_intensity', 0.1);
    parser.addParameter('nb_repetitions', 10);
    parser.addParameter('dmd_width', 1024); % px
    parser.addParameter('dmd_height', 768); % px
    parser.addParameter('dmd_frame_rate', 60.0); % Hz
    parser.addParameter('dmd_inversed_polarity', false);
    parser.addParameter('output_foldername', pwd);
    % % Parse input arguments.
    parser.parse(input_args{:});
    % % Retrieve values of input parameters.
    args = parser.Results;
    
    % Seed the random number generator.
    rng(args.seed);
    
    % Define the nodes of the gird for moving bars.
    % % Define node columns x-coordinates.
    x = (1:args.nb_node_columns);
    x = x(:) - 0.5 * (args.nb_node_columns + 1.0);
    x = x(:) * args.node_spacing;
    args.x_node_columns = x;
    % % Define node rows y-coordinate if necessary.
    y = (1:args.nb_node_rows);
    y = y(:) - 0.5 * (args.nb_node_rows + 1.0);
    y = y(:) * args.node_spacing;
    args.y_node_rows = y;
    % % Define nodes
    [X, Y] = meshgrid(args.x_node_columns, args.y_node_rows);
    
    % Plot first control_figure.
    % % ...
    x_min_60 = -7.5 * 60.0; % µm
    y_min_60 = -7.5 * 60.0; % µm
    width_60 = +15.0 * 60.0; % µm
    height_60 = +15.0 * 60.0; % µm
    x_min_30 = -7.5 * 30.0; % µm
    y_min_30 = -7.5 * 30.0; % µm
    width_30 = +15.0 * 30.0; % µm
    height_30 = +15.0 * 30.0; % µm
    % % ...
    control_1_filename = [mname, '_control_1.svg'];
    control_1_pathname = fullfile(args.output_foldername, control_1_filename);
    figure('visible', 'off');
    hold('on');
    rectangle('Position', [x_min_60, y_min_60, width_60, height_60]);
    rectangle('Position', [x_min_30, y_min_30, width_30, height_30]);
    scatter(X(:), Y(:), 50, '.');
    hold('off');
    xlim(([1.0-0.5, args.dmd_width+0.5] - (args.dmd_width + 1) / 2) * args.pixel_size);
    ylim(([1.0-0.5, args.dmd_height+0.5] - (args.dmd_height + 1) / 2) * args.pixel_size);
    title('Trace nodes over 256MEA60 and 256MEA30');
    xlabel('x (µm)');
    ylabel('y (µm)');
    saveas(gcf, control_1_pathname);
    close();
    
    % Define each trace (i.e. bar coordinates for each time step).
    args.nb_h_traces = length(args.y_node_rows);
    args.nb_v_traces = length(args.x_node_columns);
    args.nb_fd_traces = args.nb_h_traces + args.nb_v_traces - 1;
    args.nb_sd_traces = args.nb_h_traces + args.nb_v_traces - 1;
    % % Define each horizontal trace.
    h_traces = cell(args.nb_h_traces, 1);
    for h_trace_id = 1:args.nb_h_traces
        x_min = min(args.x_node_columns);
        x_max = max(args.x_node_columns);
        y_ref = args.y_node_rows(h_trace_id);
        off = args.trace_pad_length;
        x_min = x_min - off;
        x_max = x_max + off;
        d = norm([x_max, y_ref] - [x_min, y_ref]); % µm
        t = d / args.bar_speed; % sec
        n = ceil(t * args.dmd_frame_rate);
        t = n / args.dmd_frame_rate; % sec
        d = t * args.bar_speed; % µm
        x_min = (x_min + x_max - d) / 2;
        x_max = (x_min + x_max + d) / 2;
        x = linspace(x_min, x_max, n);
        y = repmat(y_ref, 1, n);
        h_traces{h_trace_id} = [x(:), y(:)];
    end
    % % Define each first diagonal trace.
    fd_traces = cell(args.nb_fd_traces, 1);
    for fd_trace_id = 1:args.nb_fd_traces
        x_min = args.x_node_columns(max(1, args.nb_v_traces - fd_trace_id + 1));
        x_max = args.x_node_columns(min(args.nb_v_traces, args.nb_fd_traces - fd_trace_id + 1));
        y_min = args.y_node_rows(max(1, fd_trace_id - args.nb_h_traces + 1));
        y_max = args.y_node_rows(min(args.nb_h_traces, fd_trace_id));
        off = args.trace_pad_length / sqrt(2);
        x_min = x_min - off;
        x_max = x_max + off;
        y_min = y_min - off;
        y_max = y_max + off;
        d = norm([x_max, y_max] - [x_min, y_min]); % µm
        t = d / args.bar_speed; % sec
        n = ceil(t * args.dmd_frame_rate);
        t = n / args.dmd_frame_rate; % sec
        d = t * args.bar_speed; % µm
        x_min = (x_min + x_max - d / sqrt(2)) / 2;
        x_max = (x_min + x_max + d / sqrt(2)) / 2;
        y_min = (y_min + y_max - d / sqrt(2)) / 2;
        y_max = (y_min + y_max + d / sqrt(2)) / 2;
        x = linspace(x_min, x_max, n);
        y = linspace(y_min, y_max, n);
        fd_traces{fd_trace_id} = [x(:), y(:)];
    end
    % % Define each vertical trace.
    v_traces = cell(args.nb_v_traces, 1);
    for v_trace_id = 1:args.nb_v_traces
        x_ref = args.x_node_columns(v_trace_id);
        y_min = min(args.y_node_rows);
        y_max = max(args.y_node_rows);
        off = args.trace_pad_length;
        y_min = y_min - off;
        y_max = y_max + off;
        d = norm([x_ref, y_max] - [x_ref, y_min]); % µm
        t = d / args.bar_speed; % sec
        n = ceil(t * args.dmd_frame_rate);
        t = n / args.dmd_frame_rate; % sec
        d = t * args.bar_speed; % µm
        y_min = (y_min + y_max - d) / 2;
        y_max = (y_min + y_max + d) / 2;
        x = repmat(x_ref, 1, n);
        y = linspace(y_min, y_max, n);
        v_traces{v_trace_id} = [x(:), y(:)];
    end
    % % Define each second diagonal trace.
    sd_traces = cell(args.nb_sd_traces, 1);
    for sd_trace_id = 1:args.nb_sd_traces
        x_min = args.x_node_columns(min(args.nb_v_traces, sd_trace_id));
        x_max = args.x_node_columns(max(1, sd_trace_id - args.nb_h_traces + 1));
        y_min = args.y_node_rows(max(1, sd_trace_id - args.nb_v_traces + 1));
        y_max = args.y_node_rows(min(args.nb_h_traces, sd_trace_id));
        off = args.trace_pad_length / sqrt(2);
        x_min = x_min + off;
        x_max = x_max - off;
        y_min = y_min - off;
        y_max = y_max + off;
        d = norm([x_max, y_max] - [x_min, y_min]); % µm
        t = d / args.bar_speed; % sec
        n = ceil(t * args.dmd_frame_rate);
        t = n / args.dmd_frame_rate; % sec
        d = t * args.bar_speed; % µm
        x_min = (x_min + x_max + d / sqrt(2)) / 2;
        x_max = (x_min + x_max - d / sqrt(2)) / 2;
        y_min = (y_min + y_max - d / sqrt(2)) / 2;
        y_max = (y_min + y_max + d / sqrt(2)) / 2;
        x = linspace(x_min, x_max, n);
        y = linspace(y_min, y_max, n);
        sd_traces{sd_trace_id} = [x(:), y(:)];
    end
    % TODO define each trace (i.e. bar coordinates for each time step).
    
    % TODO plot second control figure.
    control_2_filename = [mname, '_control_2.svg'];
    control_2_pathname = fullfile(args.output_foldername, control_2_filename);
    figure('visible', 'off');
    hold('on');
    % % Plot horizontal traces.
    for k = 1:args.nb_h_traces
        plot(h_traces{k}(:, 1), h_traces{k}(:, 2), 'b.-');
    end
    % % Plot first diagonal traces.
    for k = 1:args.nb_fd_traces
        plot(fd_traces{k}(:, 1), fd_traces{k}(:, 2), 'g.-');
    end
    % % Plot vertical traces.
    for k = 1:args.nb_v_traces
        plot(v_traces{k}(:, 1), v_traces{k}(:, 2), 'r.-');
    end
    % % Plot second diagonal traces.
    for k = 1:args.nb_sd_traces
        plot(sd_traces{k}(:, 1), sd_traces{k}(:, 2), 'k.-');
    end
    hold('off');
    xlim(([1-0.5, 1024+0.5] - (1024 + 1) / 2) * args.pixel_size);
    ylim(([1-0.5, 768+0.5] - (768 + 1) / 2) * args.pixel_size);
    xlabel('x (µm)');
    ylabel('y (µm)');
    saveas(gcf, control_2_pathname);
    
    % Define all the traces.
    normal_traces = [
        h_traces;
        fd_traces;
        v_traces;
        sd_traces
    ];
    inversed_traces = cellfun(@flipud, normal_traces, 'UniformOutput', false);
    traces = [
        normal_traces;
        inversed_traces;
    ];
    
    % Save all the traces.
    args.nb_traces = length(traces);
    trace_foldername = fullfile(args.output_foldername, mname);
    if ~exist(trace_foldername, 'dir')
        mkdir(trace_foldername)
    end
    for k = 1:args.nb_traces
        % % Open trace file.
        trace_filename = ['stimulus_', num2str(k), '.csv'];
        trace_pathname = fullfile(trace_foldername, trace_filename);
        trace_fid = fopen(trace_pathname, 'w');
        % % Write trace file header.
        trace_header = 'frameId;x;y';
        fprintf(trace_fid, '%s\r\n', trace_header);
        % % Close trace file.
        fclose(trace_fid);
        % % Write trace file data.
        trace = [(1:length(traces{k}))', traces{k}];
        dlmwrite(trace_pathname, trace, '-append', 'delimiter', ';', 'newline', 'pc');
    end
    
    % Define the schedule of the traces.
    schedule = zeros(args.nb_traces, args.nb_repetitions);
    for repetition_id = 1:args.nb_repetitions
        schedule(:, repetition_id) = randperm(args.nb_traces);
    end
    schedule = schedule(:);
    nb_trials = length(schedule);
    nb_images = cellfun(@(trace) size(trace, 1), traces);
    nb_frames = nb_images(schedule);
    start_frame_ids = cumsum([1; nb_frames(1:end-1)]);
    end_frame_ids = cumsum(nb_frames);
    orientation_ids = [
        zeros(args.nb_h_traces, 1);
        ones(args.nb_fd_traces, 1);
        repmat(2, args.nb_v_traces, 1);
        repmat(3, args.nb_sd_traces, 1);
        repmat(4, args.nb_h_traces, 1);
        repmat(5, args.nb_fd_traces, 1);
        repmat(6, args.nb_v_traces, 1);
        repmat(7, args.nb_sd_traces, 1);
    ];
    all_orientation_ids = orientation_ids(schedule);
    repetition_ids = repmat((1:args.nb_repetitions), args.nb_traces, 1);
    repetition_ids = repetition_ids(:);
    trials = [(1:nb_trials)', schedule, start_frame_ids, end_frame_ids, all_orientation_ids, repetition_ids];
    
    nb_total_images = sum(nb_images); % total number of images (i.e. unique frames)
    nb_total_frames = sum(nb_frames); % total number of frames
    
    % Write trials files.
    % % Open trials file.
    trials_filename = [mname, '_trials.csv'];
    trials_pathname = fullfile(args.output_foldername, trials_filename);
    trials_fid = fopen(trials_pathname, 'w');
    % % Write trials file header.
    trials_header = 'trialId;stimulusId;startFrameId;endFrameId;orientationId;repetitionId';
    fprintf(trials_fid, '%s\r\n', trials_header);
    % % Close trials file.
    fclose(trials_fid);
    % % Write trials file data.
    dlmwrite(trials_pathname, trials, '-append', 'delimiter', ';', 'newline', 'pc');
    
    % Write .bin file.
    % % Open .bin file.
    args.bin_filename = [mname, '_', num2str(args.dmd_frame_rate), 'hz.bin'];
    bin_pathname = fullfile(args.output_foldername, args.bin_filename);
    permission = 'w'; % writing mode, discard existing contents
    machine_format = 'l'; % IEEE floating point with little-endian byte ordering
    bin_fid = fopen(bin_pathname, permission, machine_format);
    % % Write .bin file header.
    nb_bits = 8; % number of bits
    bin_header = [args.dmd_width, args.dmd_height, nb_total_images, nb_bits];
    fwrite(bin_fid, bin_header, 'int16');
    % % Write .bin file images.
    for trace_id = 1:args.nb_traces
        trace = traces{trace_id};
        for image_id = 1:nb_images(trace_id)
            x = trace(image_id, 1);
            y = trace(image_id, 2);
            a = (orientation_ids(trace_id) / 8) * (2 * pi);
            bin_image = moving_bars_generate_image(args.dmd_width, args.dmd_height, x, y, a, args.bar_width, args.bar_length, args.pixel_size);
            bin_image = args.background_intensity + (256.0 - args.background_intensity) * bin_image;
            bin_image = uint8(bin_image - 0.5);
            bin_image = bin_image';
            if args.dmd_inversed_polarity
                fwrite(bin_fid, 255 - bin_image(:), 'uint8');
            else
                fwrite(bin_fid, bin_image(:), 'uint8');
            end
        end
    end
    % % Close .bin file.
    fclose(bin_fid);
    
    % Write .vec file.
    % % Open .vec file.
    args.vec_filename = [mname, '_', num2str(args.dmd_frame_rate), 'hz.vec'];
    vec_pathname = fullfile(args.output_foldername, args.vec_filename);
    permission = 'w'; % writing mode, discard existing contents
    machine_format = 'l'; % IEEE floating point with little-endian byte ordering
    vec_fid = fopen(vec_pathname, permission, machine_format);
    % % Write .vec file header.
    vec_header = [0, nb_total_frames, 0, 0, 0];
    fprintf(vec_fid, '%g %g %g %g %g\n', vec_header);
    % % Write .vec file image indices.
    start_image_ids = cumsum([1; nb_images(1:end-1)]);
    end_image_ids = cumsum(nb_images);
    for trial_id = 1:nb_trials
        stim_id = schedule(trial_id);
        for frame_id = start_image_ids(stim_id):end_image_ids(stim_id)
            image_id = frame_id - 1;
            vec_frame = [0, image_id, 0, 0, 0];
            fprintf(vec_fid, '%g %g %g %g %g\n', vec_frame);
        end
    end
    % % Close .vec file.
    fclose(vec_fid);
    
    % Write .mat file.
    args.mat_filename = [mname, '_parameters.mat'];
    mat_pathname = fullfile(args.output_foldername, args.mat_filename);
    save(mat_pathname, '-struct', 'args');
    
end

