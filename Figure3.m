clear all;
addpath('./myfunctions');


%% FIGURE 3 - Impact of release frequency on PDFS stability with low (default) predator mortality
disp("Computing and plotting Figure 3 - impact of release frequency on PDFS stability with low predator mortality");

%% Parameters
% Load default parameters
run("parameters.m");

% Time horizon
tspan=0:1:T;

%% Initialise figure
fig3=figure(3);
clf(fig3);
set(fig3,'Units','centimeters','Position',[0,0,30,10],'PaperUnits','centimeters','PaperSize',[30,10]*1.05);


%% FIGURE 3 (a) - stability regions

rapid_test = true;
% Specific parameter ranges
if rapid_test               % Very reduced ranges for TEST purposes
   disp("Very reduced a and Lambda_P ranges for TEST purposes");
   va=0.5*[0.1:0.5:3];      % vector of spore consumption rates 
   vLambda_P=100*[5:5:40];  % vector of yearly release quantities of predator P
else                        % Values used in Figure 3:
   va=0.5*[0.1:0.1:3];      % vector of spore consumption rates
   vLambda_P=100*[5:0.1:40];% vector of yearly release quantities of predator P
end

%% Stability of m-cPDFS for vm release frequencies & va x vLambda_P pairs
Mono_Mat=zeros(2,2); % monodromy matrix initialisation
Fond_Mat=zeros(2,2); % fondamental matrix initialisation
% stability is the PDFS stability index matrix, whose values for each (a,Lambda_P) pair are:
% -1 if unstable, +1 if stable
stability=-ones(length(va),length(vLambda_P)); % initialisation to unstable PDFS
for j=1:length(va)             % loop on spore consumption rates va
    a=va(j);                   % consumption rate
    for k=1:length(vLambda_P)  % loop on yearly release quantities of P vLambda_P
        Lambda_P=vLambda_P(k); % yearly release quantity
        for l=[0,1]            % loop on initial conditions
            y0=[1,l]-[l,0];    % initial conditions [1,0] and [0,1]
            % Integration of the linearised sub-system [I,U]:
            [tt,yy]=ode15s(@linearised_dynamics,tspan,y0);
            % Fondamental matrix:
            Fond_Mat(:,l+1)=yy(end,:)';
        end
        % Monodromy matrix:
        Mono_Mat= mtimes(diag([phi_I phi_U]),Fond_Mat);
        % Spectral radius of the monodromy matrix:
        rho=max(abs(eig(Mono_Mat)));
        % PDFS stability index being initialised to -1 standing for unstable,
        % update only to 1 if stable:
        if rho<1
           stability(j,k)=1;
        end
    end
end

%% Plot region delimitation, i.e. 0-level contour line of PDFS stability index matrix:
subplot(1,2,1)
hold on;
contour(vLambda_P, va, stability, [0 0], 'b', 'linewidth',3);

% Plot point
Lambda_P=3000;
a=0.6;
plot(Lambda_P, a, 'b.', 'LineWidth', 2, 'MarkerSize', 25);

% Subplot numbering, labels and legend
grid on; box on;
ax=gca;
ax.FontSize = 12;
text(ax.XLim(2)*0.92,ax.YLim(2)*1.04,'\bf (a)','Color','k','FontSize',12,'Interpreter','latex');
xlabel("\bf Yearly predator release~$\Lambda _P$",'Interpreter','latex');
ylabel("\bf Consumption rate~$a$",'Interpreter','latex');
text(1500,0.5,'\bf Unstable','Color','k','FontSize',12,'Interpreter','latex');
text(2300,1.1,'\bf Stable','Color','k','FontSize',12,'Interpreter','latex');
legend('$\mathcal{R}_c=1$', '$(3000, 0.6)$','Interpreter','latex')


%% FIGURE 3 (b) - stability regions

% Specific parameters
Lambda_P=3000;
a=1;

rapid_test = true;
% Specific parameter ranges
if rapid_test               % Very reduced ranges for TEST purposes
   disp("Very reduced gamma and nu ranges for TEST purposes");
   vgamma=2*[0.01:0.5:3];   % vector of sporulation rates 
   vnu=0.09*[0.1:0.2:1.3];  % vector of deposition rates
else                        % Values used in Figure 3:
   vgamma=2*[0.01:0.01:3];  % vector of sporulation rates 
   vnu=0.09*[0.1:0.005:1.3];% vector of deposition rates
end

%% Stability of PDFS and cPDFS for vgamma x vnu pairs
Mono_Mat=zeros(2,2); % monodromy matrix initialisation
Fond_Mat=zeros(2,2); % fondamental matrix initialisation
% stability[c] is the [c]PDFS stability index matrix, whose values for each (a,Lambda_P) pair are:
% -1 if unstable, +1 if stable
stability=-ones(length(vgamma),length(vnu));  % initialisation to unstable PDFS
stabilityc=-ones(length(vgamma),length(vnu)); % initialisation to unstable cPDFS
for j=1:length(vgamma)  % loop on sporulation rates vgamma
    gama=vgamma(j);     % sporulation rate
    for k=1:length(vnu) % loop on deposition rates vnu
        nu=vnu(k);      % deposition rate
        % Reproduction number determines PDFS stability:
        R=R0();
        if R<1
           stability(j,k)=1;
        end
        % cPDFS stability
        for l=[0,1]            % loop on initial conditions
            y0=[1,l]-[l,0];    % initial conditions [1,0] and [0,1]
            % Integration of the linearised sub-system [I,U]:
            [tt,yy]=ode15s(@linearised_dynamics,tspan,y0);
            % Fondamental matrix:
            Fond_Mat(:,l+1)=yy(end,:)';
        end
        % Monodromy matrix:
        Mono_Mat= mtimes(diag([phi_I phi_U]),Fond_Mat);
        % Spectral radius of the monodromy matrix:
        rho=max(abs(eig(Mono_Mat)));
        % PDFS stability index being initialised to -1 standing for unstable,
        % update only to 1 if stable:
        if rho<1
           stabilityc(j,k)=1;
        end
    end
end

%% Plot region delimitation, i.e. 0-level contour line of [c]PDFS stability index matrix:
subplot(1,2,2)
hold on; grid on; box on;
contour(vnu, vgamma, stability,  [0 0], 'k', 'linewidth',3); % PDFS
contour(vnu, vgamma, stabilityc, [0 0], 'b', 'linewidth',3); % cPDFS

% Plot point
plot(0.09, 2, '.', 'Color', [0 0.2 0.5] , 'LineWidth', 2, 'MarkerSize', 25);

% Subplot numbering, labels and legend
ax=gca;
ax.FontSize = 12;
text(ax.XLim(2)*0.92,ax.YLim(2)*1.04,'\bf (b)','Color','k','FontSize',12,'Interpreter','latex');
xlabel("\bf Deposition rate~$\nu$",'Interpreter','latex');
ylabel("\bf Sporulation rate~$\gamma$",'Interpreter','latex')
text(0.04,1,'\bf Stable','Color','k','FontSize',12,'Interpreter','latex')
text(0.04,4,'\bf Unstable','Color','k','FontSize',12,'Interpreter','latex')
legend('$\mathcal{R}_c=1$', '$\mathcal{R}=1$','$(0.09, 2)$','Interpreter','latex')


%% Save figure
saveas(fig3,'./Figure3.pdf');

 
