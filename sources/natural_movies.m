function [ ] = natural_movies( input_args )
%NATURAL_MOVIES Summary of this function goes here
%   Detailed explanation goes here
    
    [mfoldername, mname, ~] = fileparts(mfilename('fullpath'));
    
    % 1. Parse input parameters.
    % % Define input parser.
    parser = inputParser;
    parser.addParameter('input_foldername', fullfile('..', 'data'));
    parser.addParameter('input_filename', 'clips_stephanie_palmer.avi');
    parser.addParameter('video_duration', 0.0); % sec
    parser.addParameter('nb_movies', 1);
    parser.addParameter('seed', 42); % TODO remove if useless.
    parser.addParameter('background_intensity', 0.1);
    parser.addParameter('initial_adaptation_duration', 20.0); % sec
    parser.addParameter('trial_adaptation_duration', 10.0); % sec
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
    
    % TODO remove if useless.
    % Seed the random number generator.
    rng(args.seed);
    
    % ...
    input_pathname = fullfile(args.input_foldername, args.input_filename);
    disp(['Video filename: ', input_pathname]);
    video = VideoReader(input_pathname);
    if args.video_duration == 0.0
        args.video_duration = video.Duration; % sec
    end
    args.video_width = video.Width; % px
    args.video_height = video.Height; % px
    args.video_frame_rate = video.FrameRate; % Hz
    args.video_nb_bits = video.BitsPerPixel;
    args.video_nb_frames = args.video_duration * args.video_frame_rate;
    disp(['Video duration: ', num2str(args.video_duration), ' sec']);
    disp(['Number of video frames: ', num2str(args.video_nb_frames)]);
    
    % ...
    assert(mod(args.video_nb_frames, args.nb_movies) == 0);
    nb_frames = (args.video_nb_frames / args.nb_movies) * ones(args.nb_movies, 1);
    nb_frames_init = ceil(args.initial_adaptation_duration * args.dmd_frame_rate);
    nb_frames_trial = ceil(args.trial_adaptation_duration * args.dmd_frame_rate);
    args.nb_images = args.video_nb_frames + 1;
    args.nb_frames = nb_frames_init + args.nb_repetitions * (nb_frames_trial + args.video_nb_frames);
    
    % ...
    if args.background_intensity < 1.0
        args.background_intensity = 256.0 * args.background_intensity;
    else
        args.background_intensity = 1.0 * args.background_intensity;
    end
    disp(['Background intensity: ', num2str(args.background_intensity)]);
    
    % Save .bin file.
    args.bin_filename = [mname, '_', num2str(args.dmd_frame_rate),'hz.bin'];
    bin_pathname = fullfile(args.output_foldername, args.bin_filename);
    if ~exist(bin_pathname, 'file')
        % % Open .bin file.
        permission = 'w'; % writing mode, discard existing contents
        machine_format = 'l'; % IEEE floating point with little-endian byte ordering
        bin_fid = fopen(bin_pathname, permission, machine_format);
        % % Write .bin file header.
        nb_bits = 8; % number of bits
        bin_header = [args.dmd_width, args.dmd_height, args.nb_images, nb_bits];
        fwrite(bin_fid, bin_header, 'int16');
        % % Write .bin file images.
        % % % Adaptation image.
        bin_image = uint8(args.background_intensity - 0.5) * ones(args.dmd_height, args.dmd_width, 'uint8');
        bin_image = bin_image';
        if args.dmd_inversed_polarity
            fwrite(bin_fid, 255 - bin_image(:), 'uint8');
        else
            fwrite(bin_fid, bin_image(:), 'uint8');
        end
        % % % Movie images.
        i_start = floor((1 + args.dmd_height) / 2) - floor(args.video_height / 2);
        i_end = i_start + args.video_height - 1;
        j_start = floor((1 + args.dmd_width) / 2) - floor(args.video_width / 2);
        j_end = j_start + args.video_width - 1;
        for i = 1:args.video_nb_frames
            bin_image = uint8(args.background_intensity - 0.5) * ones(args.dmd_height, args.dmd_width, 'uint8');
            frame = readFrame(video);
            bin_image(i_start:i_end, j_start:j_end) = frame;
            bin_image = bin_image';
            if args.dmd_inversed_polarity
                fwrite(bin_fid, 255 - bin_image(:), 'uint8');
            else
                fwrite(bin_fid, bin_image(:), 'uint8');
            end
        end
        % % Close .bin file.
        fclose(bin_fid);
    end
    
    % Save .vec file.
    args.vec_filename = [mname, '_', num2str(args.dmd_frame_rate),'hz.vec'];
    vec_pathname = fullfile(args.output_foldername, args.vec_filename);
    if ~exist(vec_pathname, 'file')
        % % Open .vec file.
        permission = 'w'; % writing mode, discard existing contents
        machine_format = 'l'; % IEEE floating point with little-endian byte ordering
        vec_fid = fopen(vec_pathname, permission, machine_format);
        % % Write .vec file header.
        vec_header = [0, args.nb_frames, 0, 0, 0];
        fprintf(vec_fid, '%g %g %g %g %g\n', vec_header);
        % % Write .vec file image indices.
        for frame_id = 1:nb_frames_init
            image_id = 0;
            vec_frame = [0, image_id, 0, 0, 0];
            fprintf(vec_fid, '%g %g %g %g %g\n', vec_frame);
        end
        for repetition_id = 1:args.nb_repetitions
            for frame_if = 1:nb_frames_trial
                image_id = 0;
                vec_frame = [0, image_id, 0, 0, 0];
                fprintf(vec_fid, '%g %g %g %g %g\n', vec_frame);    
            end
            for frame_id = 1:args.video_nb_frames
                image_id = frame_id;
                vec_frame = [0, image_id, 0, 0, 0];
                fprintf(vec_fid, '%g %g %g %g %g\n', vec_frame);    
            end
        end
        % % Close .vec file.
        fclose(vec_fid);
    end
    
    % Define .csv file data.
    args.nb_trials = args.nb_movies * args.nb_repetitions;
    trial_ids = (1:args.nb_trials)';
    start_frame_ids = nb_frames_init + nb_frames_trial + 1 + ((1:args.nb_trials)' - 1) * (args.video_nb_frames + nb_frames_trial); % TODO correct temporary solution.
    end_frame_ids = nb_frames_init + (1:args.nb_trials)' * (args.video_nb_frames + nb_frames_trial); % TODO correct temporary solution.
    movie_ids = ones(args.nb_trials, 1); % TODO correct temporary solution.
    trials = [trial_ids, start_frame_ids, end_frame_ids, movie_ids];
    
    % Write .csv file.
    % % Open .csv file.
    csv_filename = [mname, '_trials.csv'];
    csv_pathname = fullfile(args.output_foldername, csv_filename);
    csv_fid = fopen(csv_pathname, 'w');
    % % Write .csv file header.
    csv_header = 'trialId;startFrameId;endFrameId;movieId';
    fprintf(csv_fid, '%s\r\n', csv_header);
    % % Close .csv file.
    fclose(csv_fid);
    % % Write .csv file data.
    dlmwrite(csv_pathname, trials, '-append', 'delimiter', ';', 'newline', 'pc');
    
    % Save .mat file.
    mat_filename = fullfile(args.output_foldername, [mname, '_parameters.mat']);
    save(mat_filename, '-struct', 'args');
    
end

