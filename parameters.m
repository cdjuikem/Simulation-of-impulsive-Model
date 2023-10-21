%% Default parameter values

% Global parameters
global T Lambda delta_S theta Mu mu_B mu_F omega nu gama mu_U delta_I d a e K mu_P Lambda_P m phi_S phi_I phi_U phi_P

% Parameter values
% Rainy season
T        = 245;          % Rainy season duration (days)
% Coffee leaves
Lambda   = 8;           % Recruitment rate of S
delta_S  = 0.08;        % Flower production rate by S
theta    = 0.0055;      % Maturation rate of F
Mu       = 0.0034;      % Mortality rate of leaves
mu_B     = 2.4*10^(-3); % Mortality rate of B
mu_F     = 2.4*10^(-3); % Mortality rate of F
% Rust
omega    = 0.045;       % Germination rate
nu       = 0.09;        % Deposition rate
gama     = 2;           % Sporulation rate by I
mu_U     = 0.015;       % Mortality rate of spores U
delta_I  = 0.04;        % Flower production rate by I
d        = 0.056;       % Mortality rate of I due to CLR
% Predator
a        = 0.6;         % Spore consumption rate by P
e        = 0.7;         % Biomass transformation rate
K        = 100000;      % Saturation constant of P
mu_P     = 0.003;       % Mortality rate of P
Lambda_P = 3000;        % Yearly release quantity of P
m        = 1;           % Yearly number of releases
      
% Dry season
phi_S   = 0.7;          % Survival proportion of S
phi_I   = 0.3;          % Survival proportion of I
phi_U   = 0.1;          % Survival proportion of U
phi_P   = 0.3;          % Survival proportion of P
