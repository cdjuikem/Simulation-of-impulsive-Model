clear all;
addpath('./myfunctions');


%% FIGURE 5 - Predator dynamics for cPDFS with low (default) and high predator mortality
disp("Computing and plotting Figure 5 -  Predator dynamics for cPDFS with low and high predator mortality");

%% Parameters
% Load default parameters
run("parameters.m");

% Release frequencies
vm=[1, 2, 5, 50];           % vector of release frequencies
vcolor=[0,0,1;0.6,0.4,0.3;1,0,1;0.301,0.745,0.933];% map of plot line colours (number of rows = length vm)

% Time horizon
kmax=3; % number of years simulated% Initial condition without spores

%% Initialise figure
fig5=figure(5);
clf(fig5);
set(gcf,'units','centimeters','Position',[0,0,30,15],'paperunits','centimeters','papersize',[30,15]*1.05);


%% FIGURE 5 (a) Low predator mortality

for i=1:length(vm)
    m=vm(i);

    %% Initial condition
    P0=P0_mcPDFS();
    Initial=[500,0,0,P0,0,0]; % [S(0),I(0),U(0),P(0),F(0),B(0)]

    %% System integration
    [t,u]=impulsive_solution(Initial,kmax);

    %% Plot P dynamics
    subplot(2,1,1)
    hold on; grid on; box on;
    plot(t/T,u(:,4), 'Color',colormap(vcolor(i, :)), 'linewidth',2);
    %xlabel('\bf Time (years)','Interpreter','latex');
    xticks([1:kmax]);
    title('\bf Predator $\bf P^T_m(t)$ with low mortality $\bf \mu_P=0.003/day$','Interpreter','latex')
end 


%% FIGURE 5 (b) High predator mortality
mu_P=0.1;

for i=1:length(vm)
    m=vm(i);

    %% Initial condition
    P0=P0_mcPDFS();
    Initial=[500,0,0,P0,0,0]; % [S(0),I(0),U(0),P(0),F(0),B(0)]

    %% System integration
    [t,u]=impulsive_solution(Initial,kmax);

    %% Plot P dynamics
    subplot(2,1,2)
    hold on; grid on; box on;
    plot(t/T,u(:,4), 'Color',colormap(vcolor(i, :)), 'linewidth',2);
    xlabel('\bf Time (years)','Interpreter','latex');
    xticks([1:kmax]);
    title('\bf Predator $\bf P^T_m(t)$ with high mortality $\bf \mu_P=0.1/day$','Interpreter','latex');
end 


%% Subplot numbering
subplot(2,1,1)
ax=gca;
ax.FontSize = 12;
text(ax.XLim(2)*0.92,ax.YLim(2)*1.04,'\bf (a)','Color','k','FontSize',12,'Interpreter','latex');
subplot(2,1,2)
ax=gca;
ax.FontSize = 12;
text(ax.XLim(2)*0.92,ax.YLim(2)*1.04,'\bf (b)','Color','k','FontSize',12,'Interpreter','latex');


%% Save figure
saveas(fig5,'./Figure5.pdf');


