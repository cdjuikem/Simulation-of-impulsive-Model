function P0=P0_mcPDFS()
% Compute the predator initial condition P0 for the m-CPDFS according to equation (33)
  
% Global parameters
  global T Lambda delta_S theta Mu mu_B mu_F omega nu gama mu_U delta_I d a e K mu_P Lambda_P m phi_S phi_I phi_U phi_P
  
M1=phi_P*(exp(-mu_P*T/m)-exp(-mu_P*T))+1-exp(-mu_P*T/m); % numerator of (33)
M2=(1-phi_P*exp(-mu_P*T))*(1-exp(-mu_P*T/m));            % denominator of (33)
P0plus=(Lambda_P/m)*(M1/M2);                             % equation (33)
P0=P0plus-(Lambda_P/m);

end
