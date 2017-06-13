function [ lum ] = euler_ampl_chirp( input_args )
%EULER_AMPL_CHIRP Summary of this function goes here
%   Detailed explanation goes here
    
    % Parse input parameters.
    % % Define default value for each input parameter.
    d_ante = 1.0; % sec
    d = 8.0; % sec
    d_post = 2.0; % sec
    frame_rate = 60.0; % Hz
    nb_periods = 16;
    % % Define input parser.
    parser = inputParser;
    parser.addParameter('d_ante', d_ante);
    parser.addParameter('d', d);
    parser.addParameter('d_post', d_post);
    parser.addParameter('frame_rate', frame_rate);
    parser.addParameter('nb_periods', nb_periods);
    % % Parse input arguments.
    parser.parse(input_args{:});
    % % Retrieve values of input parameters.
    args = parser.Results;
    
    
    d_total = args.d_ante + args.d + args.d_post;
    n_total = ceil(d_total * args.frame_rate);
    
    t = ((1:n_total) - 1) / args.frame_rate - args.d_ante;
    f = @(t) (t < 0.0) .* 0.5 + ((0.0 <= t) & (t < d)) .* (0.5 * (1.0 + (t / d) .* sin(2.0 * pi * args.nb_periods .* (t / d)))) + (d <= t) .* 0.5;
    lum = f(t(:));
    
end

