This repository contains Matlab code to generate DMD stimuli. They are used
with [ViALUX](https://www.vialux.de/en/index.html) products such as the
[V-Modules](https://www.vialux.de/en/v-modules.html).


## DMD Stimuli

The stimuli are:
- [`natural_full_field`](sources/natural_full_field.m) which is a natural
  full-field stimulus, its generation is based on a natural image
- [`euler_full_field`](sources/euler_full_field/euler_full_field.m) which
  is a synthetic full-field stimulus, it contains a step, a chirp in
  frequency and a chirp in amplitude
- [`flash_full_field`](source/flash_full_field.m) which is a synthetic
  stimulus, it contains flashes of light and is useful to adjust the
  maximal luminance of the light source
- [`moving_bars`](sources/moving_bars/moving_bars.m) which is a synthetic
  stimulus, bars cover a 2-D grid of nodes/points in eight directions
- [`rectangle_central_field`](source/rectangle_central_field.m) which is a synthetic
  stimulus, containing a flash in form of a rectangle of custom size and luminance 
- [`binary_time_series`](source/binary_time_series.m) which is a synthetic
  stimulus, mapping an input binary sequence into a DMD stimulus of custom luminance for 0 and 1


## Launch the code

1. Launch MATLAB.
2. Set `dmd-stimuli` as your current folder.
3. Add all folders and subfolders to path.
4. Run `initialize` in the command window.
5. Run `example_euler_full_field` in the command window.
6. Look at the outputs in the `outputs` folder.

