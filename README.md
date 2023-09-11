# FamiliarityTaxis

BAWrapper runs the experiments twice; once along the route in one direction, and once in reverse.

Bilateral_Analysis is the main function running the experiment along one route direction.
First, the images are processed to make the sky white. (WhiteSky)
Then they get cut to particular fields of view, in case of bilateral views. (split)
After this, the images get rotated and the rIDFs generated, while applying gray-scale.
