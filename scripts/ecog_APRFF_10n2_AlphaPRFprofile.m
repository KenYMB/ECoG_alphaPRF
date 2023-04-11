%% Spatial profile of alpha pRF

% 20230109
% 20230215 - update for makeFigure

% close all; clear all;
SetDefault('plotallrois',false);

hF = gobjects(0);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%
%% V1–V3
% bb_center = 20.4;  % pt
% alpha_sd1 = [16.64, 23.84]; % pt from bb_center
% 
% x0    = mean(bb_center+alpha_sd1.*[-1 1])./bb_center; % normalized alpha_center
% sigma = mean(alpha_sd1)./bb_center; % normalized alpha_sigma

%%% Parameters
x = linspace(-3,5,50000);
x0 = 1.176470588;               % center location of alpha pRF (1 for broadband)
sigma = 0.992156862745098;      % size of alpha pRF (considering true ciecle)
ss = 0.2;                       % surround suppression

%% log-x with fitting
modelfun = @(parms,x) parms(1).^(x-parms(2)) + (1-parms(1).^(x0-parms(2)));

x1 = x0 + [-1,1].*sigma;        % 1sd location at alpha pRF
[~,x1idx] = min(abs(x - x1'),[],2);
y= normpdf(x,x0,sigma)./normpdf(0,0,sigma)*-(1+ss)+ss;
parms0   = [2.2,x0];
parmslim = [0,0; inf,inf];
% parmslim = [0,x0; inf,x0];

parms = lsqcurvefit(modelfun,parms0,[1,x1],[x0,x1],parmslim(1,:),parmslim(2,:));
nx = modelfun(parms,x);

%-- Plot
hF(end+1) = figure('Color', 'w');
if plotallrois
set(gcf, 'Position', round(get(gcf,'Position').*[1/3,1,1,1]));
end
colororder(circshift(colororder(groot),-1,1));
hold on; box off;

plot(nx,y,'LineWidth',5);

[~,x1log] = min(abs(x-x1'),[],2);
x1log = nx(x1log');
hl=plot(x1log,mean(y(x1idx)).*ones(size(x1)),'k-.','LineWidth',2);
hl(2)=plot([1 1],ylim,'k:', 'LineWidth', 2); 
hl(3)=plot(xlim,[0 0],'k--','LineWidth',2); uistack(hl,'bottom');

text(1.1, -0.2, 'Center\newlinelocation', 'FontSize', 20)
text(max(x1)+0.1, mean(y(x1idx)), '1 s.d.', 'FontSize', 20)

xlim([0 5]); ylim(hl(2).YData);
box off;
set(gca, 'FontSize', 20); 
ylabel('pRF gain', 'FontSize',24)
xlabel('Eccentricity (normalized)', 'FontSize',24)

if plotallrois

%% %%%%%%%%%%%%%%%%%%%%%%%%%%
%% Dorsolateral
% bb_center = 20.4;  % pt
% alpha_sd1 = [48.31, 79.21]; % pt from bb_center
% 
% x0    = mean(bb_center+alpha_sd1.*[-1 1])./bb_center; % normalized alpha_center
% sigma = mean(alpha_sd1)./bb_center; % normalized alpha_sigma

%%% Parameters
x = linspace(-9,9,50000);
x0 = 1.176470588;               % center location of alpha pRF (1 for broadband)
sigma = 3.125490196078431;      % size of alpha pRF (considering true ciecle)
ss = 0.2;                       % surround suppression

%% log-x with fitting
modelfun = @(parms,x) parms(1).^(x-parms(2)) + (1-parms(1).^(x0-parms(2)));

x1 = x0 + [-1,1].*sigma;        % 1sd location at alpha pRF
[~,x1idx] = min(abs(x - x1'),[],2);
y= normpdf(x,x0,sigma)./normpdf(0,0,sigma)*-(1+ss)+ss;
parms0   = [2.2,x0];
parmslim = [0,0; inf,inf];
% parmslim = [0,x0; inf,x0];

parms = lsqcurvefit(modelfun,parms0,[1,x1],[x0,x1],parmslim(1,:),parmslim(2,:));
nx = modelfun(parms,x);

%-- Plot
hF(end+1) = figure('Color', 'w');
set(gcf, 'Position', round(get(gcf,'Position').*[0,1,1,1]) + (get(hF(end-1),'Position')*[1,0,1,0]').*[1,0,0,0]);
colororder(circshift(colororder(groot),-1,1));
hold on; box off;

plot(nx,y,'LineWidth',5);

[~,x1log] = min(abs(x-x1'),[],2);
x1log = nx(x1log');
hl=plot(x1log,mean(y(x1idx)).*ones(size(x1)),'k-.','LineWidth',2);
hl(2)=plot([1 1],ylim,'k:', 'LineWidth', 2); 
hl(3)=plot(xlim,[0 0],'k--','LineWidth',2); uistack(hl,'bottom');

text(1.1, -0.2, 'Center\newlinelocation', 'FontSize', 20)
text(max(x1)+1.0, mean(y(x1idx)), '1 s.d.', 'FontSize', 20)

xlim([-1 20]); ylim(hl(2).YData);
box off;
set(gca, 'FontSize', 20); 
ylabel('pRF gain', 'FontSize',24)
xlabel('Eccentricity (normalized)', 'FontSize',24)

end
