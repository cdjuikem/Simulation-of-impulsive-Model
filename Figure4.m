clear all;
addpath("./myfunctions");


%% FIGURE 4 - [S,I,U,B] dynamics with and without predators
disp("Computing and plotting Figure 4 - solution with and without predators");

%% Parameters
% Load default parameters
run("parameters.m");
% Time horizon
kmax=8; % number of years simulated

%% Initial condition
Initial=[500,0,3000,0,0,0]; % [S(0),I(0),U(0),P(0),F(0),B(0)]

%% Initialise figure
fig4=figure(4);
clf(fig4);
set(fig4,'units','centimeters','Position',[0,0,30,20],'paperunits','centimeters','papersize',[30,20]*1.05);


%% FIGURE 4 - with predators P
disp("Computing and plotting Figure 4 -  solution with predators");

%% Basic reproduction number - with P
R0_withP=R0();
disp('Reproduction number R='),disp(R0_withP);

%% System integration - with P
[t,u]=impulsive_solution(Initial,kmax);

%% Plot system dynamics - with P
color='b'; % solid blue lines
run("SIUP_plots.m");


%% FIGURE 4 - without predators
disp("Computing and plotting Figure 4 - solution without predators");

% Specific parameter(s) without predator
Lambda_P=0;

%% Basic reproduction number - without P
R0_withoutP=R0();
disp('Reproduction number R='),disp(R0_withoutP);

%% System integration - endemic
[t,u]=impulsive_solution(Initial,kmax);

%% Plot system dynamics - without P
color='k'; % solid black lines
run("SIUP_plots.m");

% Subplot numbering & legend
color='k';
subplot(2,2,1)
ax=gca;
text(ax.XLim(2)*0.92,ax.YLim(2)*1.04,'\bf (a)','Color',color,'FontSize',12,'Interpreter','latex');
subplot(2,2,2)
ax=gca;
text(ax.XLim(2)*0.92,ax.YLim(2)*1.04,'\bf (b)','Color',color,'FontSize',12,'Interpreter','latex');
subplot(2,2,3)
ax=gca;
text(ax.XLim(2)*0.92,ax.YLim(2)*1.04,'\bf (c)','Color',color,'FontSize',12,'Interpreter','latex');
subplot(2,2,4)
legend('\bf With P','\bf Without P','Interpreter','latex','Location','northeast');
ax=gca;
text(ax.XLim(2)*0.92,ax.YLim(2)*1.04,'\bf (d)','Color',color,'FontSize',12,'Interpreter','latex');


%% Save figure
saveas(fig4,'./Figure4.pdf');
