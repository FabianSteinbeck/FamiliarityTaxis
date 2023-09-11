function[sims] = rmf_split(query_img, ref_img, side, overlap, blind)

% Rotational Matching Function.
% Rotates a query image and compares it with one or more reference images
% :param query_img:
% :param ref_imgs:
% :return:

d_range = [0, 359];
d_step = 3;
degrees = [0:d_step:d_range(end)];

total_search_angles = round((d_range(2) - d_range(1)) / d_step);
sims = zeros(1, total_search_angles);

for i = 1:length(degrees)
    % rotated query image panorama
    rqimg = rotation(degrees(i), query_img);
    % rotated query image one side
    if side == 'l'
        [l,r] = split(rqimg,overlap,blind);
        rqimgs = l;
    elseif side == 'r'
        [l,r] = split(rqimg,overlap,blind);
        rqimgs = r;
    end
    sims(i) = corCoef(rqimgs, ref_img);    
end

end