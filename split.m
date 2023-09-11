function[left,right] = split(im, overlap, blind)

%% Taking a panoramic image, a left and a right image are cut with an overlap and blindspot in mind

nocs = size(im,2); % number of columns

if blind > 0
    nocspd = round(nocs/360); %nocs per degree
    bps = blind*nocspd; %blind pixels
    bpsps = bps/2; % bps per side
    im = im(:, bpsps:nocs - bpsps);
end

nocs = size(im,2); % number of columns
sp = round(nocs/2); % split point

if overlap > 0 && overlap < 360
    nocspd = round(nocs/(360-blind)); %nocs per degree
    pol = overlap * nocspd; %pixel overlap
    l_split = sp + round(pol/2);
    r_split = sp - round(pol/2);
    left = im(:, 1:l_split);
    right = im(:, r_split:nocs);
elseif overlap == 360
    left = im;
    right = im;
elseif overlap < 0
    nocspd = round(nocs/(360-blind)); %nocs per degree
    pol = overlap * nocspd; %pixel overlap
    l_split = sp + round(pol/2);
    r_split = sp - round(pol/2);
    left = im(:, 1:l_split);
    right = im(:, r_split:nocs);
else
    left = im(:, 1:sp);
    right = im(:, sp+1:nocs);
end

end