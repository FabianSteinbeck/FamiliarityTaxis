function[im_w] = whiteSky(im)

% Turns grey sky white by identifying the most pixels which have the same values (the sky is completely of one value)
% :param query_img:
% :param ref_imgs:
% :return:

max_val = 255; %white

[c,bL] = imhist(im(:,:,1)); % most values are sky

max_im = find(c == max(c));
sky_val = bL(max_im);

im_w = im;
im_w(im == sky_val) = max_val; % replace sky with white

end