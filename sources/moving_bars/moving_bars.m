function [ ] = moving_bars( input_args )
%MOVING_BARS Summary of this function goes here
%   Detailed explanation goes here
    
    [~, mname, ~] = fileparts(mfilename('fullpath'));
    
    % 1. Parse input parameters.
    % % Define default value for each input parameter.
    pixel_size = 2.3; % µm
    background_intensity = 0.1;
    nb_repetitions = 5 * 30;
    dmd_width = 1024; % px
    dmd_height = 768; % px
    dmd_frame_rate = 60.0; % Hz
    dmd_inversed_polarity = false;
    output_foldername = pwd;
    % % Define input parser.
    parser = inputParser;
    parser.addParameter('pixel_size', pixel_size);
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
    
    % Plot first control_figure.
    A = zeros(args.dmd_height, args.dmd_width, 3, 'uint8');
    x_min = -7.5 * 60.0; % µm
    x_max = +7.5 * 60.0; % µm
    y_min = -7.5 * 60.0; % µm
    y_max = +7.5 * 60.0; % µm
    j_min = ceil(0.5 * args.dmd_width + x_min / args.pixel_size);
    j_max = ceil(0.5 * args.dmd_width + x_max / args.pixel_size);
    i_min = ceil(0.5 * args.dmd_height - y_max / args.pixel_size);
    i_max = ceil(0.5 * args.dmd_height - y_min / args.pixel_size);
    A(i_min:i_max, j_min:j_max, :) = 127;
    x_min = -7.5 * 30.0; % µm
    x_max = +7.5 * 30.0; % µm
    y_min = -7.5 * 30.0; % µm
    y_max = +7.5 * 30.0; % µm
    j_min = ceil(0.5 * args.dmd_width + x_min / args.pixel_size);
    j_max = ceil(0.5 * args.dmd_width + x_max / args.pixel_size);
    i_min = ceil(0.5 * args.dmd_height - y_max / args.pixel_size);
    i_max = ceil(0.5 * args.dmd_height - y_min / args.pixel_size);
    A(i_min:i_max, j_min:j_max, :) = 255;
    
    x = ceil((-2:+2) * 300.0 / args.pixel_size + 0.5 * args.dmd_width);
    y = ceil((-2:+2) * 300.0 / args.pixel_size + 0.5 * args.dmd_height);
    [X, Y] = meshgrid(x, y);
    
    control_1_filename = [mname, '_control_1.svg'];
    control_1_pathname = fullfile(args.output_foldername, control_1_filename);
    figure('visible', 'off');
    hold('on');
    image(A);
    scatter(X(:), Y(:));
    hold('off');
    saveas(gcf, control_1_pathname);
    close();
    
end

