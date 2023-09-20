# FamiliarityTaxis

BAWrapper runs the experiments twice; once along the route in one direction, and once in reverse.

Bilateral_Analysis is the main function running the experiment along one route direction. (_N for the standard resolution, _Compressed for the various, downsampled images)

First, the images are processed to make the sky white. (WhiteSky)

Then they get cut to particular fields of view, in case of bilateral views. (split)

After this, the images get rotated (rotation) and the rIDFs generated while applying gray-scale. (rmf_split, cor_coef).

Once all rIDFs are generated, the differences between left and right eyes are calculated (the off-route images normalised). From those differences, the correct steering values are identified (left of heading values must be negative, right of heading must be positive). These then are summed up in different ways for plotting.

The Aliasing is analysed in the AliasAnalysis, which is very similar to the Bilateral_Analysis with a few extra steps towards the end, calculating the difference between the identified location along the route vs the actual position along the route. 
