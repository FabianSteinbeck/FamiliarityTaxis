function[l,r] = offset(img, offset)

% Rotational Matching Function.
% Rotates a query image and compares it with one or more reference images
% :param query_img:
% :param ref_imgs:
% :return:

d_step_l = -offset; % left eye rotates right 
d_step_r = offset; % right eye rotates left

l = rotation(d_step_l,img);
r = rotation(d_step_r,img);

end