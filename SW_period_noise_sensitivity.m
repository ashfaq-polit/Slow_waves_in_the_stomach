close all;
clear all;

% % 
% load('SW_period_267_42_3point5CaIP340X_change_both.mat'); a = SW_period
% load('SW_period_267_42_5point3CaIP340X_change_both.mat'); b = SW_period
% load('SW_period_267_42_7CaIP340X_change_both.mat');       c = SW_period
% load('SW_period_267_42_14CaIP340X_change_both.mat');      d = SW_period

load('Total_lag_267_42_7CaIP30_half_both_gjdev2510_ip3dev2510.mat');  
a = lag1(end-6:end);b = lag2(end-6:end);c = lag3(end-6:end);

x = 1:1:3;
y = [mean(a) mean(b) mean(c)];
err = [std(a) std(b) std(c)];
errorbar(x,y,err,'vertical','s','Markersize',2,'MarkerFaceColor','r',...
    'MarkerEdgeColor','red','LineWidth',2)

ax = gca; box off;
ax.LineWidth = 3;
ax.FontSize = 22;
ax.FontName = 'arial';
% ylim([20.5 21.5])
xlim([0 4])
xlabel('% Variability');ylabel('SW(period)')

All_period = [a; b; c]'
[p,tbl,stats] = anova1(All_period);
[h2,p2]= ttest2(b,c,'Vartype','unequal')

% legend('0.35 nS','0.53 nS','0.70 nS','1.4 nS')
% legend('4.0/sec','6.0/sec','8.0/sec','16.0/sec')
xlabel('Cycle #');
% ylabel('SW_period (sec)')
box off; hold off

figure
violin(All_period)
ylim([0 105])
%% Multi-variable analysis

number = 41
gj_array = zeros(number,3);
for i=1:number
    gj_array(i,1) = .70 + (-1+2*rand(1,1))*.2*.7
    gj_array(i,2) = .70 + (-1+2*rand(1,1))*.5*.7
    gj_array(i,3) = .70 + (-1+2*rand(1,1))*1.0*.7
end

number = 41
gj_array = zeros(number,3);
for i=1:number
    gj_array(i,1) =  (-1+2*rand(1,1))*.2
    gj_array(i,2) =  (-1+2*rand(1,1))*.5
    gj_array(i,3) =  (-1+2*rand(1,1))*1.0
end
    
%// Define x values
x = (1:41).';
xMat = repmat(x, 1, 3); %// For plot3

%// Define y values
y = 1:1:3;
yMat = repmat(y, numel(x), 1); %//For plot3

zMat = [gj_array(:,1) gj_array(:,2) gj_array(:,3)]; %// For plot3

plot3(yMat, xMat,zMat, '.-','MarkerSize',15); %// Make all traces blue
line([3,3],[0,41],[0,0],'Color','k','Linewidth',3) 
ax = gca; box off;
ax.LineWidth = 3;
ax.FontSize = 22;
ax.FontName = 'arial';
ylim([1 41])
set(gcf, 'Renderer', 'painters');
set(gcf,'Position',[100 100 700 500])



%% Bar chart attempt

%https://stackoverflow.com/questions/24987216/stacking-multiple-2d-plots-into-a-single-3d-plot-in-matlab
%https://www.mathworks.com/matlabcentral/answers/373286-plot-multiple-variables-in-different-colors-with-scatter3

number = 41
gj_array = zeros(number,3);
for i=1:number
    gj_array(i,1) = .70 + (-1+2*rand(1,1))*.2*.7
    gj_array(i,2) = .70 + (-1+2*rand(1,1))*.5*.7
    gj_array(i,3) = .70 + (-1+2*rand(1,1))*1.0*.7
end

%// Define x values
x = (1:41).';
xMat1 = repmat(x, 1, 1); %// For plot3

%// Define y values
y = 1:1:3;
yMat1 = repmat(y(1), numel(x), 1); %//For plot3

zMat1 = [gj_array(:,1)]; %// For plot3


xMat2 = repmat(x, 1, 1); %// For plot3
yMat2 = repmat(y(2), numel(x), 1); %//For plot3
zMat2 = [gj_array(:,2)]; %// For plot3


xMat3 = repmat(x, 1, 1); %// For plot3
yMat3 = repmat(y(3), numel(x), 1); %//For plot3
zMat3 = [gj_array(:,3)]; %// For plot3

figure
plot3(yMat1, xMat1, zMat1, 'b.-','MarkerSize',15); %// Make all traces blue
hold on
plot3(yMat2, xMat2, zMat2,'m.-','MarkerSize',15); %
hold on
plot3(yMat3, xMat3, zMat3, 'g.-','MarkerSize',15); %

ax = gca; box off;
ax.LineWidth = 3;
ax.FontSize = 22;
ax.FontName = 'arial';
ylim([1 41])
set(gcf, 'Renderer', 'painters');
set(gcf,'Position',[100 100 700 500])

