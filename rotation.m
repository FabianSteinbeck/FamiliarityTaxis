function[rotIm] = rotation(angle,im)
% Converts the degrees into columns and rotates the image around the vertical axis.
% Positive degrees turn the image anti-clockwise (left)
% and negative degrees rotate the image clockwise (right)
% angle: number of degrees the agent will rotate its view
% im: An array that we want to shift.
% rotIm: Returns the rotated image.
%%
nocs = size(im,2);
nocspd = nocs/360;
cols_to_shift = round(angle * nocspd);
rotIm = circshift(im, cols_to_shift, 2);

end