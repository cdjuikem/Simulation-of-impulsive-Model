clear all;
addpath('./myfunctions');


%% FIGURE 6 - Impact of release frequency on m-cPDFS stability with low (default) predator mortality
disp("Computing and plotting Figure 6 - impact of release frequency on m-cPDFS stability with low predator mortality");

%% Parameters
% Load default parameters
run("parameters.m");

%% Initialise figure
fig6=figure(6);
clf(fig6);
set(fig6,'Units','centimeters','Position',[0,0,30,10],'PaperUnits','centimeters','PaperSize',[30,10]*1.05);


%% FIGURE 6 (a) - stability regions

rapid_test = true;
% Specific parameter ranges
if rapid_test               % Very reduced ranges for TEST purposes
   disp("Very reduced a and Lambda_P ranges for TEST purposes");
   va=0.5*[0.1:0.5:3];      % vector of spore consumption rates 
   vLambda_P=100*[5:5:40];  % vector of yearly release quantities of predator P
else                        % Values used in Figure 6:
   va=0.5*[0.1:0.1:3];      % vector of spore consumption rates
   vLambda_P=100*[5:0.1:40];% vector of yearly release quantities of predator P
end
vm=[1, 5];                  % vector of release frequencies
vcolor=['b', 'm'];          % vector of plot line colours (same length as vm)

% Time horizon
tspan=0:1:T;

subplot(1,2,1)
hold on;

%% Stability of m-cPDFS for vm release frequencies & va x vLambda_P pairs
for i=1:length(vm) % loop on release frequencies
    m=vm(i);
    color=vcolor(i);
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
    contour(vLambda_P, va, stabilityN, [0 0], color, 'linewidth',3);
end

% Plot point for dynamics subplot
Lambda_P=3000;
a=0.6;
plot(Lambda_P, a, 'c.', 'LineWidth', 2, 'MarkerSize', 25);

% Subplot labels, legend and numbering
grid on; box on;
ax=gca;
ax.FontSize = 12;
text(ax.XLim(2)*0.92,ax.YLim(2)*1.04,'\bf (a)','Color','k','FontSize',12,'Interpreter','latex');
xlabel("\bf Yearly predator release~$\Lambda _P$",'Interpreter','latex');
ylabel("\bf Consumption rate~$a$",'Interpreter','latex');
text(1500,0.5,'\bf Unstable','Color','k','Interpreter','latex');
text(2300,1.2,'\bf Stable','Color','k','Interpreter','latex');
legend('$m=1$','$m=5$','$(3000, 0.6)$','Interpreter','latex');


%% FIGURE 6 (b) - I dynamics

% Basic reproduction number
R0=R0();
disp('Reproduction number R='),disp(R0);

% Time horizon
kmax=12; % number of years simulated

% Initial condition
Initial=[500,0,3000,0,0,0]; % [S(0),I(0),U(0),P(0),F(0),B(0)]

subplot(1,2,2)
hold on;

for i=1:length(vm) % loop on release frequencies
    m = vm(i);
    color= vcolor(i);

    %% System integration
    [t,u]=impulsive_solution(Initial,kmax);

    %% Plot I dynamics
    plot(t/T,u(:,2),color,'linewidth',2);
end  

% Subplot numbering and labels
grid on; box on;
ax=gca;
text(ax.XLim(2)*0.92,ax.YLim(2)*1.04,'\bf (b)','Color','k','FontSize',12,'Interpreter','latex');
ax.FontSize = 12;
xlabel('\bf Time (years)','Interpreter','latex');
xticks([1:kmax]);
title('\bf Infected leaves I','Interpreter','latex'); 


%% Save figure
saveas(fig6,'./Figure6.pdf');

 
