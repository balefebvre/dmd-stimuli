function show_trace(lum, num_samples, fs)
% Function to display a part of stimulus in its final shape.

rng default
N = numel(lum);
start = randi([1 N-num_samples-1]); % pick a random part of the stimulus
x = lum(start:start+num_samples);

figure, title('Stimulus example')
plot((floor(1:numel(x))/fs), x, 'k', 'LineWidth', 0.6), box off
xlabel('Time (s)'), ylabel('Luminosity (0-255)')

end

