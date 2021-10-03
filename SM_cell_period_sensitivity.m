close all;
clear all;

load('Period_267_42_3point55point3710point5CaIP340X_half_both.mat'); 
a = period_last1; b = period_last2; c = period_last3; d = period_last4;
load('Period_267_42_14203050CaIP340X_half_both.mat');
e = period_last1; f = period_last2; g = period_last3; h = period_last4;

load('Period_267_42_7CaIP320X30X40X60X_half_both.mat'); 
a = period_last1; b = period_last2; c = period_last3; d = period_last4;
load('Period_267_42_7CaIP380X100X150X200X_half_both.mat');
e = period_last1; f = period_last2; g = period_last3; h = period_last4;

load('Period_267_42_7Cazeropoint05point1point2IP340X_half_both'); 
a = period_last1; b = period_last2; c = period_last3; d = period_last4;
load('Period_267_42_7Capoint51510IP340X_half_both.mat');
e = period_last1; f = period_last2; g = period_last3; h = period_last4;

%% Calculations
SW_periodu1 = a; SW_periodu2 = b; SW_periodu3 = c; SW_periodu4 = d;
SW_periodu5 = e; SW_periodu6 = f; SW_periodu7 = g; SW_periodu8 = h;


mean_SW = [mean(SW_periodu1(length(SW_periodu1)-6:length(SW_periodu1))) mean(SW_periodu2(length(SW_periodu2)-6:length(SW_periodu2)))...
    mean(SW_periodu3(length(SW_periodu3)-6:length(SW_periodu3))) mean(SW_periodu4(length(SW_periodu4)-6:length(SW_periodu4)))...
    mean(SW_periodu5(length(SW_periodu5)-6:length(SW_periodu5))) mean(SW_periodu6(length(SW_periodu6)-6:length(SW_periodu6)))...
    mean(SW_periodu7(length(SW_periodu7)-6:length(SW_periodu7))) mean(SW_periodu8(length(SW_periodu8)-6:length(SW_periodu8)))];
%     mean(SW_periodu9(length(SW_periodu9)-6:length(SW_periodu9)))];
std_SW = [std(SW_periodu1(length(SW_periodu1)-6:length(SW_periodu1))) std(SW_periodu2(length(SW_periodu2)-6:length(SW_periodu2)))...
    std(SW_periodu3(length(SW_periodu3)-6:length(SW_periodu3))) std(SW_periodu4(length(SW_periodu4)-6:length(SW_periodu4)))...
    std(SW_periodu5(length(SW_periodu5)-6:length(SW_periodu5))) std(SW_periodu6(length(SW_periodu6)-6:length(SW_periodu6)))...
    std(SW_periodu7(length(SW_periodu7)-6:length(SW_periodu7))) std(SW_periodu8(length(SW_periodu8)-6:length(SW_periodu8)))]
%     std(SW_periodu9(length(SW_periodu9)-6:length(SW_periodu9)))];
figure
plot(SW_periodu1(length(SW_periodu1)-6:length(SW_periodu1)),'r.-', 'MarkerSize', 25,'LineWidth',2); 
hold on; plot(SW_periodu2(length(SW_periodu2)-6:length(SW_periodu2)),'g.-', 'MarkerSize', 25, 'LineWidth', 2); 
hold on; plot(SW_periodu3(length(SW_periodu3)-6:length(SW_periodu3)),'b.-', 'MarkerSize', 25, 'LineWidth', 2); 
hold on; plot(SW_periodu4(length(SW_periodu4)-6:length(SW_periodu4)),'c.-', 'MarkerSize', 25, 'LineWidth', 2); 
hold on; plot(SW_periodu5(length(SW_periodu5)-6:length(SW_periodu5)),'m.-', 'MarkerSize', 25, 'LineWidth', 2); 
hold on; plot(SW_periodu6(length(SW_periodu6)-6:length(SW_periodu6)),'y.-', 'MarkerSize', 25, 'LineWidth', 2); 
hold on; plot(SW_periodu7(length(SW_periodu7)-6:length(SW_periodu7)),'k.-', 'MarkerSize', 25, 'LineWidth', 2); 
line([1,7],[0,0],'Color','k')
ylim([17 20]) 

figure
% plot(0.35,SW_periodu1(length(SW_periodu1)-6:length(SW_periodu1)),'k.','Markersize',25)
% hold on; 
plot(0.53,SW_periodu2(length(SW_periodu2)-6:length(SW_periodu2)),'k.','Markersize',25)
hold on; plot(0.7,SW_periodu3(length(SW_periodu3)-6:length(SW_periodu3)),'k.','Markersize',25)
hold on; plot(1.05,SW_periodu4(length(SW_periodu4)-6:length(SW_periodu4)),'k.','Markersize',25)
hold on; plot(1.4,SW_periodu5(length(SW_periodu5)-6:length(SW_periodu5)),'k.','Markersize',25)
hold on; plot(2.0,SW_periodu6(length(SW_periodu6)-6:length(SW_periodu6)),'k.','Markersize',25)
hold on; plot(3.0,SW_periodu7(length(SW_periodu7)-6:length(SW_periodu7)),'k.','Markersize',25)
hold on; plot(5.0,SW_periodu8(length(SW_periodu8)-6:length(SW_periodu8)),'k.','Markersize',25)

figure
plot(0,SW_periodu1(length(SW_periodu1)-6:length(SW_periodu1)),'k.','Markersize',25)
hold on; 
plot(0.05,SW_periodu2(length(SW_periodu2)-6:length(SW_periodu2)),'k.','Markersize',25)
hold on; plot(0.1,SW_periodu3(length(SW_periodu3)-6:length(SW_periodu3)),'k.','Markersize',25)
hold on; plot(0.2,SW_periodu4(length(SW_periodu4)-6:length(SW_periodu4)),'k.','Markersize',25)
hold on; plot(0.5,SW_periodu5(length(SW_periodu5)-6:length(SW_periodu5)),'k.','Markersize',25)
hold on; plot(1.0,SW_periodu6(length(SW_periodu6)-6:length(SW_periodu6)),'k.','Markersize',25)
hold on; plot(5.0,SW_periodu7(length(SW_periodu7)-6:length(SW_periodu7)),'k.','Markersize',25)
hold on; plot(10.0,SW_periodu8(length(SW_periodu8)-6:length(SW_periodu8)),'k.','Markersize',25)

ax = gca; box off;
ax.LineWidth = 3;
ax.FontSize = 22;
ax.FontName = 'arial';
xlim([0 10.0]);
ylim([17 18])
xlabel('% Variability');ylabel('SM Cell Period')
hold on

% x=[8.0 12.0 16.0 20.0 30.0 40.0]'
% x = [0.53 0.7 1.05 1.4 2.0 3.0 5.0]'
x=[0 0.05 0.1 0.2 0.5 1.0 5.0 10.0]';
y = [mean(SW_periodu1(length(SW_periodu1)-6:length(SW_periodu1)))...
    mean(SW_periodu2(length(SW_periodu2)-6:length(SW_periodu2)))...
mean(SW_periodu3(length(SW_periodu3)-6:length(SW_periodu3)))... 
mean(SW_periodu4(length(SW_periodu4)-6:length(SW_periodu4)))...
mean(SW_periodu5(length(SW_periodu5)-6:length(SW_periodu5)))...
mean(SW_periodu6(length(SW_periodu6)-6:length(SW_periodu6)))...
mean(SW_periodu7(length(SW_periodu7)-6:length(SW_periodu7)))...
mean(SW_periodu8(length(SW_periodu8)-6:length(SW_periodu8)))]';

myfittype = fittype('a1*exp(-x/beta1)+a2*exp(-x/beta2)',...
      'dependent',{'y'},'independent',{'x'},...
      'coefficients',{'a1','beta1','a2','beta2'})
myfit = fit(x,y,myfittype, 'StartPoint',[17 1 1 1])
%'StartPoint',[30 1 1 1]
coeff=coeffvalues(myfit)
plot(myfit,x,y,'k.'); 
