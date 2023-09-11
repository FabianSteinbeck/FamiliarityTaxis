function[mAmNegL,mAmPosR] = Bilateral_Analysis_N(mode)

%% Setting things up
pics = dir('img*'); % images from the antworld

% visual fields with overlap (column 1) and blindspot (column 2)
combinations = [-180,0;...
                -90,90;...
                0,180;...
                -90,0;...
                0,90;...
                90,180;...
                0,0;...
                90,90;...
                180,180;...                
                90,0;...
                180,90;...
                180,0;...
                270,90;...
                270,0;...
                360,0];

% direction of movement pull relative to agent heading
os = linspace(0,90,11);

%% Organization
% load images
im = cell(size(pics,1),1);
for i = 1:size(pics,1)
    im{i} = imread(pics(i).name);
    im{i} = whiteSky(im{i}); %make the sky white
    if mode == 1 % make them reverse
        im{i} = rotation(180,im{i});
    end
end

% dedicated cells for left, centre, right images
imResh = reshape(im,[11,size(im,1)/11]);
imL = {}; % image Left
imL = imResh(1:5,:); %flip to make row 1 being the smallest distance to centre
imC = {}; % image on Centre route
imC = imResh(6,:);
imR = {}; % image Right
imR = imResh(7:end,:);
imOC = [imL;imR]; % off-centre images

%% 360 rIDF

for i = 1:length(imC)
    rIDFP(i,:) = rotation(180,rmf(imC{i},imC{i})); % Panorama
    for ii = 1:5
        rIDFPL{ii,i} = rotation(180,rmf(imL{ii,i},imC{i})); % 360 left
        rIDFPR{ii,i} = rotation(180,rmf(imR{ii,i},imC{i})); % 360 right
        disp(strcat(num2str(i),' ',num2str(ii)))
    end
end
if mode == 1 % reverse
    save('RrIDFP','rIDFP')
    save('RrIDFPL','rIDFPL')
    save('RrIDFPR','rIDFPR')
else
    save('FrIDFP','rIDFP')
    save('FrIDFPL','rIDFPL')
    save('FrIDFPR','rIDFPR')
end

%% Bilateral rIDF
for i = 1:length(os)
    for ii = 1:size(combinations,1)
        for iii = 1:length(imC)
            % Snapshot centre split
			[memO_l,memO_r] = split(imC{iii},combinations(ii,1),combinations(ii,2));
            
            % rotate the centre-route snapshots by an offset
			[off_l,off_r] = offset(imC{iii},os(i));
            
			% rotary Image Difference Function Bilateral 
            rIDFB(i,ii,iii).l = rotation(180,rmf_split(off_l,memO_l,...
                                          'l',combinations(ii,1),combinations(ii,2)));
            rIDFB(i,ii,iii).r = rotation(180,rmf_split(off_r,memO_r,...
                                          'r',combinations(ii,1),combinations(ii,2)));
										  
			% compare off-centre query images with centre-snapshots		  
            for iv = 1:size(imOC,1)
				
				% rotate the off-route query images by an offset
                [offoff_l,offoff_r] =...
                    offset(imOC{iv,iii},os(i));
					
				% left eye	
                rIDF(i,ii,iii,iv).l = rotation(180,rmf_split(offoff_l,...
                    memO_l,'l',combinations(ii,1),combinations(ii,2))); 
					
				% right eye	
                rIDF(i,ii,iii,iv).r = rotation(180,rmf_split(offoff_r,...
                    memO_r,'r',combinations(ii,1),combinations(ii,2))); 
                disp(strcat(num2str(i),'/',num2str(length(os)),...
                    '_',num2str(ii),'/',num2str(size(combinations,1)),...
                    '_',num2str(iii),'/',num2str(length(imC)),...
                    '_',num2str(iv),'/',num2str(size(imOC,1))))
            end
        end
    end
end

if mode == 1 % reverse
    save('RrIDF','rIDF')
    save('RrIDFB','rIDFB')
else
    save('FrIDF','rIDF')
    save('FrIDFB','rIDFB')
end

%% Analysis 360

for i = 1:size(rIDFPL,1)
    for ii = 1:size(rIDFPL,2)
        maxPL(i,ii) = max(rIDFPL{i,ii});
        maxPR(i,ii) = max(rIDFPR{i,ii});
        maxPLI(i,ii) = find(rIDFPL{i,ii} == maxPL(i,ii)); %max - index
        maxPRI(i,ii) = find(rIDFPR{i,ii} == maxPR(i,ii)); %max - index
    end
end

%% Analysis Bilateral

for i = 1:size(rIDF,1) % offset
    for ii = 1:size(rIDF,2) % FoV
        for iii = 1:size(rIDF,3) %Snapshots
            
			% subtracting the right from the left rIDF for the on-route snapshots
			differenceC{i,ii,iii} = rIDFB(i,ii,iii).l -...
                                 rIDFB(i,ii,iii).r;
								 
            for iv = 1:size(rIDF,4) % Off route
                maxBL(i,ii,iii,iv) = max(rIDF(i,ii,iii,iv).l); %max val left
                maxBR(i,ii,iii,iv) = max(rIDF(i,ii,iii,iv).r); %max val right
                maxBLI(i,ii,iii,iv) = find(rIDF(i,ii,iii,iv).l...
                                        == maxBL(i,ii,iii,iv)); %max - index
                maxBRI(i,ii,iii,iv) = find(rIDF(i,ii,iii,iv).r...
                                        == maxBR(i,ii,iii,iv)); %max - index
                % 5%min-max-normalisation
                minL = min(rIDF(i,ii,iii,iv).l);
                maxL = max(rIDF(i,ii,iii,iv).l);
                l_5 = diff([minL,maxL])*0.05;
                rIDF_L = (rIDF(i,ii,iii,iv).l - (minL + l_5))/...
                            (maxL - (minL + l_5));
                minR = min(rIDF(i,ii,iii,iv).r);
                maxR = max(rIDF(i,ii,iii,iv).r);
                r_5 = diff([minR,maxR])*0.05;
                rIDF_R = (rIDF(i,ii,iii,iv).r - (minR + r_5))/...
                            (maxR - (minR + r_5));
                difference{i,ii,iii,iv} = rIDF_L - rIDF_R;
%                 difference{i,ii,iii,iv} = rIDF(i,ii,iii,iv).l -...
%                                        rIDF(i,ii,iii,iv).r;

                %leftward rotation frontal 180 deg (route following direction +- 90°)
                posL(i,ii,iii,iv) = sum(difference{i,ii,iii,iv}...
                    (round(length(difference{i,ii,iii,iv})/4):...
                    round(length(difference{i,ii,iii,iv})/2)) >=0);%leftward steering forward left quadrant
                amPosL(i,ii,iii,iv) = posL(i,ii,iii,iv)/(length(difference{i,ii,iii,iv})/4);% percent
                negL(i,ii,iii,iv) = sum(difference{i,ii,iii,iv}...
                    (round(length(difference{i,ii,iii,iv})/4):...
                    round(length(difference{i,ii,iii,iv})/2)) <0);%rightward steering forward right quadrant
                amNegL(i,ii,iii,iv) = negL(i,ii,iii,iv)/(length(difference{i,ii,iii,iv})/4);% percent

                %rightward rotation frontal 180 deg (route following direction +- 90°)
                posR(i,ii,iii,iv) = sum(difference{i,ii,iii,iv}...
                                        (round(length(difference{i,ii,iii,iv})/2)+1:...
                                        round(length(difference{i,ii,iii,iv})*3/4)) >=0); %leftward steering
                amPosR(i,ii,iii,iv) = posR(i,ii,iii,iv)/(length(difference{i,ii,iii,iv})/4);% percent
                negR(i,ii,iii,iv) = sum(difference{i,ii,iii,iv}...
                                        (round(length(difference{i,ii,iii,iv})/2)+1:...
                                        round(length(difference{i,ii,iii,iv})*3/4)) <0);%rightward steering
                amNegR(i,ii,iii,iv) = negR(i,ii,iii,iv)/(length(difference{i,ii,iii,iv})/4);% percent

                %leftward rotation frontal 360 deg (route following direction +- 180°)
                posLF(i,ii,iii,iv) = sum(difference{i,ii,iii,iv}...
                    (1:...
                    round(length(difference{i,ii,iii,iv})/2)) >=0);%leftward steering forward left quadrant
                amPosLF(i,ii,iii,iv) = posLF(i,ii,iii,iv)/(length(difference{i,ii,iii,iv})/2);% percent
                negLF(i,ii,iii,iv) = sum(difference{i,ii,iii,iv}...
                    (1:...
                    round(length(difference{i,ii,iii,iv})/2)) <0);%rightward steering forward right quadrant
                amNegLF(i,ii,iii,iv) = negLF(i,ii,iii,iv)/(length(difference{i,ii,iii,iv})/2);% percent

                %rightward rotation frontal 360 deg (route following direction +- 180°)
                posRF(i,ii,iii,iv) = sum(difference{i,ii,iii,iv}...
                                        (round(length(difference{i,ii,iii,iv})/2)+1:...
                                        round(length(difference{i,ii,iii,iv}))) >=0); %leftward steering
                amPosRF(i,ii,iii,iv) = posRF(i,ii,iii,iv)/(length(difference{i,ii,iii,iv})/2);% percent
                negRF(i,ii,iii,iv) = sum(difference{i,ii,iii,iv}...
                                        (round(length(difference{i,ii,iii,iv})/2)+1:...
                                        round(length(difference{i,ii,iii,iv}))) <0);%rightward steering
                amNegRF(i,ii,iii,iv) = negRF(i,ii,iii,iv)/(length(difference{i,ii,iii,iv})/2);% percent

            end
        end
        % 180
        mAmPosL(i,ii) = mean(amPosL(i,ii,:,:),'all');
        mAmNegL(i,ii) = mean(amNegL(i,ii,:,:),'all'); 
        mAmPosR(i,ii) = mean(amPosR(i,ii,:,:),'all');
        mAmNegR(i,ii) = mean(amNegR(i,ii,:,:),'all');
        % 360
        mAmPosLF(i,ii) = mean(amPosLF(i,ii,:,:),'all');
        mAmNegLF(i,ii) = mean(amNegLF(i,ii,:,:),'all'); 
        mAmPosRF(i,ii) = mean(amPosRF(i,ii,:,:),'all');
        mAmNegRF(i,ii) = mean(amNegRF(i,ii,:,:),'all');
    end
end

if mode == 1 % reverse
    save('FNL','mAmNegL')
    save('FPR','mAmPosR')    
else
    save('RNL','mAmNegL')
    save('RPR','mAmPosR')    
end
%% Intermediate plotting
Plt(imC,rIDFP,rIDFPL,rIDFPR,rIDF,rIDFB,difference,differenceC,...
    maxBL,maxBR,maxBLI,maxBRI,amNegL,amPosR,mAmNegL,mAmPosR,mAmNegLF,mAmPosRF,combinations,os,mode);

%% saving

%% Plotting


%%
end