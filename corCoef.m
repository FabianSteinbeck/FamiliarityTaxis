function[cor_coef] = corCoef(a, b)
%Gaussian blur, edge detection and image resize
% :param a: picture 1
% :param b: picture 2
% :return: cor_coef

a = double(a);
am = mean(a,3); % collapse rgb
b = double(b);
bm = mean(b,3);
% flatten
ar = reshape(am.',1,[]);
br = reshape(bm.',1,[]);

cor_coef = cov(ar, br) / (std(ar) * std(br));
cor_coef = cor_coef(1,2);
end