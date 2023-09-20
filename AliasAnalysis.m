dataNames = dir('*rIDF*');

r = 1;

for c = 1:size(dataNames,1)

load(dataNames(c).name)
% Analysis Bilateral

for i = 1:size(rIDF,1) % offset
    for ii = 1:size(rIDF,2) % FoV
        for iii = 1:size(rIDF,3) %Off route
            for iv = 1:size(rIDF,4) % On route
                maxBL(i,ii,iii,iv) = max(rIDF(i,ii,iii,iv).l); %max val left
                maxBR(i,ii,iii,iv) = max(rIDF(i,ii,iii,iv).r); %max val right
%                 maxBLI(i,ii,iii,iv) = find(rIDF(i,ii,iii,iv).l...
%                                         == maxBL(i,ii,iii,iv)); %max - index
%                 maxBRI(i,ii,iii,iv) = find(rIDF(i,ii,iii,iv).r...
%                                         == maxBR(i,ii,iii,iv)); %max - index
                % 5%min-max-normalisation
                minL = min(rIDF(i,ii,iii,iv).l);
                maxL = max(rIDF(i,ii,iii,iv).l);
                l_5 = diff([minL,maxL])*0.05;
                rIDF_L(iv,:) = (rIDF(i,ii,iii,iv).l - (minL + l_5))/...
                            (maxL - (minL + l_5));
                minR = min(rIDF(i,ii,iii,iv).r);
                maxR = max(rIDF(i,ii,iii,iv).r);
                r_5 = diff([minR,maxR])*0.05;
                rIDF_R(iv,:) = (rIDF(i,ii,iii,iv).r - (minR + r_5))/...
                            (maxR - (minR + r_5));
            end
            
            % on-route memory index most similar to current position on off-route
            IL = find(maxBL(i,ii,iii,:) == max(maxBL(i,ii,iii,:)));
            IR = find(maxBR(i,ii,iii,:) == max(maxBR(i,ii,iii,:)));
            
            errorL(iii) = abs(iii - IL);
            errorR(iii) = abs(iii - IR);
            
            difference{i,ii,iii} = rIDF_L(IL,:) - rIDF_R(IR,:);
%                 difference{i,ii,iii,iv} = rIDF(i,ii,iii,iv).l -...
%                                        rIDF(i,ii,iii,iv).r;
            %leftward rotation frontal 180 deg
            posL(i,ii,iii) = sum(difference{i,ii,iii}...
                (round(length(difference{i,ii,iii})/4):...
                round(length(difference{i,ii,iii})/2)) >=0);%leftward steering forward left quadrant
            amPosL(i,ii,iii) = posL(i,ii,iii)/(length(difference{i,ii,iii})/4);% percent
            negL(i,ii,iii) = sum(difference{i,ii,iii}...
                (round(length(difference{i,ii,iii})/4):...
                round(length(difference{i,ii,iii})/2)) <0);%rightward steering forward right quadrant
            amNegL(i,ii,iii) = negL(i,ii,iii)/(length(difference{i,ii,iii})/4);% percent

            %rightward rotation frontal 180 deg
            posR(i,ii,iii) = sum(difference{i,ii,iii}...
                                    (round(length(difference{i,ii,iii})/2)+1:...
                                    round(length(difference{i,ii,iii})*3/4)) >=0); %leftward steering
            amPosR(i,ii,iii) = posR(i,ii,iii)/(length(difference{i,ii,iii})/4);% percent
            negR(i,ii,iii) = sum(difference{i,ii,iii}...
                                    (round(length(difference{i,ii,iii})/2)+1:...
                                    round(length(difference{i,ii,iii})*3/4)) <0);%rightward steering
            amNegR(i,ii,iii) = negR(i,ii,iii)/(length(difference{i,ii,iii})/4);% percent

            %leftward rotation frontal 360 deg
            posLF(i,ii,iii) = sum(difference{i,ii,iii}...
                (round(length(difference{i,ii,iii})):...
                round(length(difference{i,ii,iii})/2)) >=0);%leftward steering forward left quadrant
            amPosLF(i,ii,iii) = posLF(i,ii,iii)/(length(difference{i,ii,iii})/2);% percent
            negLF(i,ii,iii) = sum(difference{i,ii,iii}...
                (round(length(difference{i,ii,iii})):...
                round(length(difference{i,ii,iii})/2)) <0);%rightward steering forward right quadrant
            amNegLF(i,ii,iii,iv) = negLF(i,ii,iii)/(length(difference{i,ii,iii})/2);% percent

            %rightward rotation frontal 360 deg
            posRF(i,ii,iii) = sum(difference{i,ii,iii}...
                                    (round(length(difference{i,ii,iii})/2)+1:...
                                    round(length(difference{i,ii,iii}))) >=0); %leftward steering
            amPosRF(i,ii,iii) = posRF(i,ii,iii)/(length(difference{i,ii,iii})/2);% percent
            negRF(i,ii,iii) = sum(difference{i,ii,iii}...
                                    (round(length(difference{i,ii,iii})/2)+1:...
                                    round(length(difference{i,ii,iii}))) <0);%rightward steering
            amNegRF(i,ii,iii) = negRF(i,ii,iii)/(length(difference{i,ii,iii})/2);% percent

        end
        
        eL(i,ii) = mean(errorL);
        eR(i,ii) = mean(errorR);
        
        % 180
        mAmNegL(i,ii) = mean(amNegL(i,ii,:),'all'); 
        mAmPosR(i,ii) = mean(amPosR(i,ii,:),'all');
        % 360
        mAmNegLF(i,ii) = mean(amNegLF(i,ii,:),'all'); 
        mAmPosRF(i,ii) = mean(amPosRF(i,ii,:),'all');
        
   end
end

%% Figures

finalSteering(:,:,c) = mean(cat(ndims(mAmNegL) + 1, mAmNegL, mAmPosR), ndims(mAmNegL) + 1);

finalAlias(:,:,c) = mean(cat(ndims(eL) + 1, eL, eR), ndims(eL) + 1);
maxFA(c) = max(finalAlias(:,:,c),[],'all');
minFA(c) = min(finalAlias(:,:,c),[],'all');

figure('units','normalized','outerposition',[0 0 1 1])
imagesc(finalAlias(:,:,c))
title('Mean Alias')
xlabel('Field of View')
xticks([1:1:15])
xticklabels({'-180/0','-90/90','0/180','-90/0','0/90','90/180','0/0',...
    '90/90','180/180','90/0','180/90','180/0','270/90','270/0','360/0'})
ylabel('Offset')
yticks([1:1:11])
yticklabels({'0','9','18','27','36','45','54','63','72','81','90'})
colorbar
colormap gray;

name = strcat(ZeroNamer('AliasRoute',c,max(size(dataNames,1))),'.png');
print(gcf,name,'-dpng','-r300');
name = strcat(ZeroNamer('AliasRoute',c,max(size(dataNames,1))));
saveas(gcf,name,'epsc')
close all

figure('units','normalized','outerposition',[0 0 1 1])
imagesc(finalSteering(:,:,c),[0.5,1])
title('Mean Steering')
xlabel('Field of View')
xticks([1:1:15])
xticklabels({'-180/0','-90/90','0/180','-90/0','0/90','90/180','0/0',...
    '90/90','180/180','90/0','180/90','180/0','270/90','270/0','360/0'})
ylabel('Offset')
yticks([1:1:11])
yticklabels({'0','9','18','27','36','45','54','63','72','81','90'})
colorbar
colormap gray;

name = strcat(ZeroNamer('SteeringRoute',c,max(size(dataNames,1))),'.png');
print(gcf,name,'-dpng','-r300');
name = strcat(ZeroNamer('SteeringRoute',c,max(size(dataNames,1))));
saveas(gcf,name,'epsc')
close all

%% saving
if mod(c,2) ~= 0 % reverse
    name = strcat(ZeroNamer('Route',r,10),'REL');
    save(name,'eL')
    name = strcat(ZeroNamer('Route',r,10),'RER');
    save(name,'eR')
    name = strcat(ZeroNamer('Route',r,10),'RNL');
    save(name,'mAmNegL')
    name = strcat(ZeroNamer('Route',r,10),'RPR');
    save(name,'mAmPosR')   
else
    name = strcat(ZeroNamer('Route',r,10),'FEL');
    save(name,'eL')
    name = strcat(ZeroNamer('Route',r,10),'FER');
    save(name,'eR')
    name = strcat(ZeroNamer('Route',r,10),'FNL');
    save(name,'mAmNegL')
    name = strcat(ZeroNamer('Route',r,10),'FPR');
    save(name,'mAmPosR')
    r = r + 1;
end

clearvars -except c dataNames r finalSteering finalAlias maxFA minFA
end

%%

fS = mean(finalSteering,3);
save('MeanSteering','fS')
for i = 1:size(finalSteering,3)
    FS75(i) = sum(finalSteering(:,:,i) > 0.75,'all');
end
save('Percent75','FS75')

figure('units','normalized','outerposition',[0 0 1 1])
imagesc(fS,[0.5,1])
title('Mean Steering')
xlabel('Field of View')
xticks([1:1:15])
xticklabels({'-180/0','-90/90','0/180','-90/0','0/90','90/180','0/0',...
    '90/90','180/180','90/0','180/90','180/0','270/90','270/0','360/0'})
ylabel('Offset')
yticks([1:1:11])
yticklabels({'0','9','18','27','36','45','54','63','72','81','90'})
colorbar
colormap gray;

name = strcat('MeanSteering','.png');
print(gcf,name,'-dpng','-r300');
name = 'MeanSteering';
saveas(gcf,name,'epsc')
close all

fA = mean(finalAlias,3);
save('MeanAlias','fA')

figure('units','normalized','outerposition',[0 0 1 1])
imagesc(fA)
title('Mean Alias')
xlabel('Field of View')
xticks([1:1:15])
xticklabels({'-180/0','-90/90','0/180','-90/0','0/90','90/180','0/0',...
    '90/90','180/180','90/0','180/90','180/0','270/90','270/0','360/0'})
ylabel('Offset')
yticks([1:1:11])
yticklabels({'0','9','18','27','36','45','54','63','72','81','90'})
colorbar
colormap gray;

name = strcat('MeanAlias','.png');
print(gcf,name,'-dpng','-r300');
name = 'MeanAlias';
saveas(gcf,name,'epsc')
close all

figure('units','normalized','outerposition',[0 0 1 1])
plot(maxFA(1:2:20),'b')
hold on
plot(minFA(1:2:20),'c')
plot(maxFA(2:2:20),'r')
plot(minFA(2:2:20),'m')

name = strcat('AliasVariation','.png');
print(gcf,name,'-dpng','-r300');
name = 'AliasVariation';
saveas(gcf,name,'epsc')
close all