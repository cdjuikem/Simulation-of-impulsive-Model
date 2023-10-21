function dyp=linearised_dynamics(t,y)
% Continuous dynamics of linearised sub-system [I,U] around the PDFS

% Global parameters
global T Lambda delta_S theta Mu mu_B mu_F omega nu gama mu_U delta_I d a e K mu_P Lambda_P m phi_S phi_I phi_U phi_P

dyp=zeros(2,1); % initialisation

dyp(1)=-(Mu+d)*y(1)+omega*nu*y(2);               % linear dynamics of infected leaves I
dyp(2)=gama*y(1)-(nu+mu_U+a*predator(t)/K)*y(2); % linear dynamics of uredospores U

end
