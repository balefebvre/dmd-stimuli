This repository contains Matlab code to generate DMD stimuli. They are used
with [ViALUX](https://www.vialux.de/en/index.html) products such as the
[V-Modules](https://www.vialux.de/en/v-modules.html).

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
- ...
