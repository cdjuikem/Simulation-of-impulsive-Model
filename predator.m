function P=predator(t)
% m-cPDFS predator value at time t according to equation (32)

% Global parameters
global T Lambda delta_S theta Mu mu_B mu_F omega nu gama mu_U delta_I d a e K mu_P Lambda_P m phi_S phi_I phi_U phi_P

M1=phi_P*(exp(-mu_P*T/m)-exp(-mu_P*T))+1-exp(-mu_P*T/m); % numerator of (33)
M2=(1-phi_P*exp(-mu_P*T))*(1-exp(-mu_P*T/m));            % denominator of (33)
P0plus=(Lambda_P/m)*(M1/M2);                             % equation (33)

j=floor(m*t/T); % number of the last release
M3=1-exp(-mu_P*j*T/m);
M4=1-exp(-mu_P*T/m);
P=(P0plus*exp(-mu_P*j*T/m)+(Lambda_P/m)*(M3/M4)).*exp(-mu_P*(t-(j*T/m))); % equation (32)

end
