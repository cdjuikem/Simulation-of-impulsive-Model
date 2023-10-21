function [tt,yp]=impulsive_solution(y0,kmax)
% Returns the solution of the impulsive system with initial condition y during kmax seasons

% Global parameters
global T Lambda delta_S theta Mu mu_B mu_F omega nu gama mu_U delta_I d a e K mu_P Lambda_P m phi_S phi_I phi_U phi_P
% Local parameter
n=100; % number of output steps between two releases (or per year without predator)

% Initialisation of variables
yp=[];
tt=[];
opts = odeset('RelTol',1e-6,'AbsTol',1e-6); % ode45 tolerance options
for k=0:kmax-1        % season loop
    for j=0:m-1       % predator release loop
      y0plus=y0+[0, 0, 0, Lambda_P/m, 0, 0]; % predator release
        dt=T/(n*m);   % output time step
        [tj,yj]=ode45(@continuous_dynamics,[(k*T+j*(T/m)):dt:(k*T+(j+1)*(T/m))],y0plus,opts); % ODE solution between releases j and j+1
        yp=[yp;yj];   % solution concatenation
        tt=[tt;tj];   % time concatenation
        y0=yp(end,:); % initial condition for next release
        if j==m-1
          y0=[phi_S, phi_I, phi_U, phi_P, 0, 0].*y0; % initial condition for next season
        end
    end
end
end
