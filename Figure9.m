clear all;
addpath('./myfunctions');
    

%% FIGURE 9 - Impact of release frequency on transient dynamics with high mortality
disp("Computing and plotting Figure 9 - impact of release frequency on transient dynamics with high predator mortality");

%% Parameters
% Load default parameters
run("parameters.m");

% Specific parameter
mu_P = 0.1; % high predator mortality

%% Initialise figure
fig9=figure(9);
clf(fig9);
set(fig9,'Units','centimeters','Position',[0,0,30,13],'PaperUnits','centimeters','PaperSize',[30,13]*1.05);


%% FIGURE 9 (a) - stability regions

rapid_test = true;
% Specific parameter ranges
if rapid_test                 % Very reduced ranges for TEST purposes
   disp("Very reduced a and Lambda_P ranges for TEST purposes");
   va=[0.1:0.5:3];            % vector of spore consumption rates 
   vLambda_P=1000*[10:10:50]; % vector of yearly release quantities of predator P
else                          % Values used in Figure 9:
   va=[0.1:0.1:3];            % vector of spore consumption rates
   vLambda_P=1000*[10:0.1:50];% vector of yearly release quantities of predator P
end
vm=[1, 2, 5, 50];             % vector of release frequencies
vcolor=[0,0,1;0.6,0.4,0.3;1,0,1;0.301,0.745,0.933];% map of plot line colours (number of rows = length vm)

% Time horizon
tspan=0:1:T;

subplot(2,2,[1,3])
hold on;
 
%% Stability of m-cPDFS for vm release frequencies & va x vLambda_P pairs
for i=1:length(vm) % loop on release frequencies
    m=vm(i);
    color=colormap(vcolor(i,:));
    Mono_MatN=zeros(2,2); % monodromy matrix initialisation
    Fond_Mat=zeros(2,2);  % fondamental matrix initialisation
    % stabilityN is the m-cPDFS stability index matrix, whose values for each (a,Lambda_P) pair are:
    % -1 if unstable, +1 if stable
    stabilityN=-ones(length(va),length(vLambda_P)); % initialisation to unstable m-cPDFS
    for j=1:length(va)             % loop on spore consumption rates va
        a=va(j);                   % consumption rate
        for k=1:length(vLambda_P)  % loop on yearly release quantities of P vLambda_P
            Lambda_P=vLambda_P(k); % yearly release quantity
            for l=[0,1]            % loop on initial conditions
                y0=[1,l]-[l,0];    % initial conditions [1,0] and [0,1]
                % Integration of the linearised sub-system [I,U]:
                [tN,yN]=ode15s(@linearised_dynamics,tspan,y0);
                % Fondamental matrix:
                Fond_Mat(:,l+1)=yN(end,:)';
            end
            % Monodromy matrix:
            Mono_MatN= mtimes(diag([phi_I phi_U]),Fond_Mat);
            % Spectral radius of the monodromy matrix:
            rho=max(abs(eig(Mono_MatN)));
            % m-cPDFS stability index being initialised to -1 standing for unstable,
            % update only to 1 if stable:
            if rho<1
                stabilityN(j,k)=1;
            end
        end
    end
    %% Plot region delimitation, i.e. 0-level contour line of m-cPDFS stability index matrix:
   contour(vLambda_P, va, stabilityN, [0 0], 'Color',color, 'linewidth',3);
end

% Plot point for dynamics subplot
Lambda_P=30000;
a=2.5;
plot(Lambda_P, a,'.','Color',[0.7 0.3 0.9],'LineWidth', 2, 'MarkerSize', 25);

% Subplot numbering, labels and legend
grid on; box on;
ax=gca;
ax.FontSize = 12;
text(ax.XLim(2)*0.92,ax.YLim(2)*1.04,'\bf (a)','Color','k','FontSize',12,'Interpreter','latex');
xlabel("\bf Yearly predator release~$\Lambda _P$",'Interpreter','latex');
ylabel("\bf Consumption rate~$a$",'Interpreter','latex');
text(18000,0.8,'\bf Unstable','Color','k','Interpreter','latex');
text(33000,2.3,'\bf Stable','Color','k','Interpreter','latex');
legend('$m=1$','$m=2$','$m=5$','$m=50$','$(30000, 2.5)$','Interpreter','latex');


%% FIGURE 9 (b) - [I,B] dynamics

% Basic reproduction number
R0=R0();
disp('Reproduction number R='),disp(R0);

% Time horizon
kmax=6; % number of years simulated

% Initial condition
Initial=[500,0,3000,0,0,0]; % [S(0),I(0),U(0),P(0),F(0),B(0)]

for i=1:length(vm) % loop on release frequencies
    m = vm(i);
    color=colormap(vcolor(i,:));

    %% ODE system integration
    [t,u]=impulsive_solution(Initial,kmax);

    % Plot I dynamics
    subplot(2,2,2)
    hold on;
    plot(t/T,u(:,2),'Color',color,'linewidth',2);

    % Plot B dynamics
    subplot(2,2,4)
    hold on;
    plot(t/T,u(:,6),'Color',color,'linewidth',2);
end  

%% Subplot numbering and labels
subplot(2,2,2)
grid on; box on;
ax=gca;
ax.FontSize = 12;
text(ax.XLim(2)*0.92,ax.YLim(2)*1.04,'\bf (b)','Color','k','FontSize',12,'Interpreter','latex');
xlabel('\bf Time (years)','Interpreter','latex')
xticks([1:kmax]);
title('\bf Infected leaves I','Interpreter','latex'); 

subplot(2,2,4)
grid on; box on;
ax=gca;
ax.FontSize = 12;
xlabel('\bf Time (years)','Interpreter','latex')
xticks([1:kmax]);
title('\bf Berries B','Interpreter','latex'); 


%% Save figure
saveas(fig9,'./Figure9.pdf');

 
