function[RMSE] = rmse(a, b)
%Image Differencing Function RMSE
% :param a: A single query image
% :param b: One or more reference images
% :return:

if size(b,3) > 0
    for ref_img = 1:size(b,3) 
        RMSE(ref_img) = sqrt(mean(ref_img - a));
    end
end

RMSE = sqrt(mean(b - a));

end