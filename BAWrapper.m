function[] = BAWrapper()

% Forward Correct Steering Response measures
[FNL,FPR] = Bilateral_Analysis_N(0);
% Reverse Correct Steering Response measures
[RNL,RPR] = Bilateral_Analysis_N(1);

%% Final plot

F = [FNL + FPR + RNL + RPR]./4;

figure('units','normalized','outerposition',[0 0 1 1])
imagesc(F > 0.75,[0.5,1])
title('Mean Correct Steering')
xlabel('Field of View')
xticks([1:1:15])
xticklabels({'-180/0','-90/90','0/180','-90/0','0/90','90/180','0/0',...
    '90/90','180/180','90/0','180/90','180/0','270/90','270/0','360/0'})
ylabel('Offset')
yticks([1:1:11])
yticklabels({'0','9','18','27','36','45','54','63','72','81','90'})
colorbar
colormap gray;

name = strcat('MeanSteeringResponse75','.png');
print(gcf,name,'-dpng','-r300');
name = strcat('MeanSteeringResponse75');
saveas(gcf,name,'epsc')
close all

%%
end