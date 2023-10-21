clear all
addpath("./myfunctions")


%% FIGURE 2 - [S,I,U,B] dynamics without predators

%% Parameters
% Load default parameters
run("parameters.m");
% No predators
Lambda_P=0;
% Time horizon
kmax=6;      % number of years simulated

%% Initial condition
Initial=[500,0,3000,0,0,0]; % [S,I,U,P,F,B]

%% Initialise figure
fig2=figure(2)
clf(fig2)
set(fig2,'units','centimeters','Position',[0,0,30,20],'paperunits','centimeters','papersize',[30,20]*1.05)


%% FIGURE 2 - endemic solution
disp("Computing and plotting Figure 2 - endemic solution without predators");

%% Basic reproduction number - endemic
R0_endemic=R0();
disp(R0_endemic)

%% System integration - endemic
[t,u]=impulsive_solution(Initial,kmax);

%% Plot system dynamics - endemic
color='k'; % solid black lines
run("SIUB_plots.m");


%% FIGURE 2 - periodic disease free solution (PDFS)
disp("Computing and plotting Figure 2 - periodic disease free solution without predators (PDFS)");

% Specific PDFS parameter(s)
gama=1.6;

%% Basic reproduction number - PDFS
R0_PDFS=R0();
disp(R0_PDFS)

%% System integration - PDFS
[t,u]=impulsive_solution(Initial,kmax);

%% Plot system dynamics - PDFS
color='-.k'; % dash-dotted black lines
run("SIUB_plots.m");

% Subplot numbering & legend
color='k';
subplot(2,2,1)
ax=gca;
text(ax.XLim(2)*0.92,ax.YLim(2)*1.04,'\bf (a)','Color',color,'FontSize',12,'Interpreter','latex')
subplot(2,2,2)
ax=gca;
text(ax.XLim(2)*0.92,ax.YLim(2)*1.04,'\bf (b)','Color',color,'FontSize',12,'Interpreter','latex')
subplot(2,2,3)
ax=gca;
text(ax.XLim(2)*0.92,ax.YLim(2)*1.04,'\bf (c)','Color',color,'FontSize',12,'Interpreter','latex')
subplot(2,2,4)
legend('$\bf \mathcal{R}>1$','$\bf \mathcal{R}<1$' ,'Interpreter','latex','Location','northwest')
ax=gca;
text(ax.XLim(2)*0.92,ax.YLim(2)*1.04,'\bf (d)','Color',color,'FontSize',12,'Interpreter','latex')


%% Save figure
saveas(fig2,'./Figure2.pdf')
