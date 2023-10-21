function dyp=continuous_dynamics(t,y)
% Continuous dynamics of system [S,I,U,P,F,B]

% Global parameters
global T Lambda delta_S theta Mu mu_B mu_F omega nu gama mu_U delta_I d a e K mu_P Lambda_P m phi_S phi_I phi_U phi_P

infec=(omega*nu*y(3))/(y(1)+y(2)); % force of infection

dyp(1)=Lambda-infec*y(1)-Mu*y(1);                     % dynamics of susceptible leaves S
dyp(2)=infec*y(1)-(Mu+d)*y(2);                        % dynamics of infected leaves I
dyp(3)=gama*y(2)-(nu+mu_U)*y(3)-a*y(3)*y(4)/(K+y(3)); % dynamics of uredospores U
dyp(4)=e*a*y(3)*y(4)/(K+y(3))-mu_P*y(4);              % dynamics of predators P
dyp(5)=delta_S*y(1)+delta_I*y(2)-(theta+mu_F)*y(5);   % dynamics of flowers F
dyp(6)=theta*y(5)-mu_B*y(6);                          % dynamics of berries B
dyp=[dyp(1);dyp(2);dyp(3);dyp(4);dyp(5);dyp(6)];      
end

