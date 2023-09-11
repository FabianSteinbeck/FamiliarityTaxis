function[sims] = rmf(query_img, ref_imgs)

% Rotational Matching Function.
% Rotates a query image and compares it with one or more reference images
% :param query_img:
% :param ref_imgs:
% :return:

d_range = [0, 359];
d_step = 3;
degrees = [0:d_step:d_range(end)];

total_search_angles = round((d_range(2) - d_range(1)) / d_step);
% sims = zeros(size(ref_imgs,3), total_search_angles);

for i = 1:length(degrees)
    % rotated query image
    rqimg = rotation(degrees(i), query_img);
    sims(i) = corCoef(rqimg, ref_imgs);    
end

end