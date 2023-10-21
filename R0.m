function R_0=R0()
% Compute R0 for the impulsive system according to equation (18)

% Global parameters
global T Lambda delta_S theta Mu mu_B mu_F omega nu gama mu_U delta_I d a e K mu_P Lambda_P m phi_S phi_I phi_U phi_P

k1=Mu+d;
k2=nu+mu_U;
alpha=k1 +k2;
beta=sqrt((k1-k2)^2+4*gama*omega*nu);
lam1=-(alpha+beta)/2;
lam2=(-alpha+beta)/2;
R_0=(phi_I/(2*beta))*((beta+k1-k2)*exp(lam1*T)+(beta-k1+k2)*exp(lam2*T));

end
