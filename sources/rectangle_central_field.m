function [ ] = rectangle_central_field( input_args )
% RECTANGLE_CENTRAL_FIELD Generate a rectangle stimulus in 
% the center of the DMD display.
    
    [mfoldername, mname, ~] = fileparts(mfilename('fullpath'));
    
    % 1. Parse input parameters.
    % % Define input parser.
    parser = inputParser;
    parser.addParameter('input_foldername', fullfile(mfoldername, '..', 'data'));
    parser.addParameter('rectangle_width', 100) % px
    parser.addParameter('rectangle_height', 50) % px
    parser.addParameter('rectangle_duration', 0.5); % sec
    parser.addParameter('post_rect_duration', 0.5); % sec
    parser.addParameter('initial_adaptation_duration', 0.0);
    parser.addParameter('trial_adaptation_duration', 1.5);
    parser.addParameter('rectangle_intensity', 1.0);
    parser.addParameter('background_intensity', 0.1);
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
          
    % Define stimulus and post-stimulus duration.
    n_rect = ceil(args.rectangle_duration * args.dmd_frame_rate);
    rectangle = ones(n_rect, 1); 
    n_post_rect = ceil(args.post_rect_duration * args.dmd_frame_rate);
    post_rectangle = zeros(n_post_rect, 1);
    
    % Composing stimulus of rectangle and post-rectangle.
    lum_stim = [
        rectangle;
        post_rectangle;
    ];
    
    % Get rectangle coordinates. 
    width = args.rectangle_width; height = args.rectangle_height;
    dmd_width = args.dmd_width; dmd_height = args.dmd_height;
    ind_x = dmd_width/2 - floor(width/2) : dmd_width/2 + ceil(width/2) - 1;
    ind_y = dmd_height/2 - floor(height/2) : dmd_height/2 + ceil(height/2) - 1;

    % Compute stimulus duration.
    n_stim = length(lum_stim);
    d_stim = n_stim / args.dmd_frame_rate;
    disp('Stimulus duration:');
    disp(['  ', num2str(d_stim), ' sec']);
    if d_stim > 60.0
        dm = d_stim / 60.0; % minutes
        ds = d_stim - dm * 60.0; % seconds
        disp(['  ', dm, ' min ', ds, ' sec']);
    end
    
    % % Rescale stimulus luminance profile.
    if args.background_intensity < 1.0
        args.background_intensity = 256.0 * args.background_intensity;
    end
    disp(['Background intensity: ', num2str(args.background_intensity)]);
    lum_stim = lum_stim * (256.0 - args.background_intensity);
    lum_stim = lum_stim + args.background_intensity;
    lum_stim = uint8(lum_stim - 0.5);
    
    if args.rectangle_intensity < 1.0
        args.rectangle_intensity = 256.0 * args.rectangle_intensity;
    end
    disp(['Rectangle intensity: ', num2str(args.rectangle_intensity)]);
     
    % Compute luminance statistics.
    disp(['Minimal luminance: ', num2str(min(lum_stim))]);
    disp(['Median luminance: ', num2str(median(lum_stim))]);
    disp(['Maximal luminance: ', num2str(max(lum_stim))]);
   
    % Generate full-field profile used at the beginning.
    n_init = ceil(args.initial_adaptation_duration * args.dmd_frame_rate);
    lum_init = args.background_intensity * ones(n_init, 1);
    lum_init = uint8(lum_init - 0.5);
    
    % Generate full-field profile used before each repetition.
    n_trial = ceil(args.trial_adaptation_duration * args.dmd_frame_rate);
    lum_trial = args.background_intensity * ones(n_trial, 1);
    lum_trial = uint8(lum_trial - 0.5);
    
    % Generate full-field profile.
    profiles = cell(1 + 2 * args.nb_repetitions, 1);
    profiles{1} = lum_init;
    for repetition_id = 1:args.nb_repetitions
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
    for i = 1:nb_images
        if mod(i,2) % background
            bin_image = zeros(args.dmd_width, args.dmd_height); 
            bin_image = args.background_intensity + (256.0 - args.background_intensity) * bin_image;
            bin_image = uint8(bin_image - 0.5);   
            if args.dmd_inversed_polarity
                fwrite(bin_fid, 255 - bin_image(:), 'uint8');
            else
                fwrite(bin_fid, bin_image(:), 'uint8');
            end
        else        % rectangle
            tmp = zeros(width, height);
            tmp = args.rectangle_intensity + (256.0 - args.rectangle_intensity) * tmp;
            bin_image(ind_x, ind_y) = uint8(tmp - 0.5);
            % for const color: bin_image(ind_x, ind_y) = 255*ones(width, height, 'uint8');
            if args.dmd_inversed_polarity
                fwrite(bin_fid, 255 - bin_image(:), 'uint8');
            else
                fwrite(bin_fid, bin_image(:), 'uint8');
            end          
        end
        % imshow(bin_image)
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
    repetitions = zeros(args.nb_repetitions, 3);
    nbs_frames = cellfun(@(lum) length(lum), profiles);
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
    
end