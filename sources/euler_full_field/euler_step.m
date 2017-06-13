function [ lum ] = euler_step( input_args )
%EULER_STEP Summary of this function goes here
%   Detailed explanation goes here
    
    % Parse input parameters.
    % % Define default value for each input parameter.
    d_ante = 1.5; % sec
    d = 3.0; % sec
    d_post = 3.0; % sec
    frame_rate = 60.0; % Hz
    % % Define input parser.
    parser = inputParser;
    parser.addParameter('d_ante', d_ante);
    parser.addParameter('d', d);
    parser.addParameter('d_post', d_post);
    parser.addParameter('frame_rate', frame_rate);
    % % Parse input arguments.
    parser.parse(input_args{:});
    % % Retrieve values of input parameters.
    args = parser.Results;
    
    
    d_total = args.d_ante + args.d + args.d_post;
    n_total = ceil(d_total * args.frame_rate);
    
    lum = zeros(n_total, 1);
    
    i_start = (1 + floor(args.d_ante * args.frame_rate)) + 1;
    i_end = (1 + floor((args.d_ante + args.d) * args.frame_rate));
    
    lum(i_start: i_end) = 1.0;
    
    return
    
end

