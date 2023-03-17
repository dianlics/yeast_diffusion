% parameters for simulation
param.tmax = 400; %maximal yeast growth duration, unit: h
param.dt = 0.1; %unit: h
param.L = 10; %maximal radius, unit: mm
param.N = 100; %grid number

% cell growth on Met
param.n_max = 10^4;
param.dn = 10^-3; %unit: mm^2/h
param.miu_max = 0.5; %unit: h^(-1)
param.K_Met = 0.6; %unit: 10^-10 mol

% Met
param.dMet = 10^-3; %unit: mm^2/h
param.vmax_Met = 6*10^-5; %unit: 10^-10 mol/hr
param.YMet_n = 10^-5; %unit: 10^-10 mol

%%%%%%%%%%%%%%%%%%%%%%%%%%%%to be determined!!!
% inhibition of cell growth from Pb2+
param.KI_Pb = 1*10^0; %unit: 10^-10 mol
param.mI_Pb = 1;

% inhibition of cell growth from S2-
param.KI_S2 = 5*10^5; %unit: 10^-10 mol
param.mI_S2 = 4;

% Met5
param.YMet5_n = 3000;
param.vmax_Met5 = 1000; %unit: h^(-1)
param.KI_Met = 30; %unit: 10^-10 mol
param.mI_Met = 4;
param.delta_Met5 = 0.01;

% S2
param.dS2 = 10^-3; %unit: mm^2/h
param.vmax_S2 = 10^-3; %unit: 10^-10 mol
param.Ksp = 1*10^-8; %unit: 10^-10 mol

% gene expression capacity parameters
param.Kphi=1.820769438301 ; %0-20
param.exp_phi=3; %1-5 randi(5)
param.phi=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% range of variable
% param.n_min = 0;
% param.Met_min = 0;
% param.Met_max = 0.134;
% param.Met5_min = 0;
% param.Met5_max = param.YMet5_n*param.n_max;
% param.SO4_min = 0;
% param.SO4_max = 1.514;
% param.S2_min = 0;
% param.S2_max = 10;
% param.Pb2_min = 0;
% param.Pb2_max = 0.00307;
% param.PbS_min = 0;
% param.PbS_max = 0.00307;

% Initial distributions
% xx=linspace(0,param.L,param.N)';
% param.n0= 10*exp(-(xx).^2/.005);
param.n0= zeros(param.N,1);
param.n0(1:10) = 100;
param.Met0 = zeros(param.N,1);
param.Met50 = param.YMet5_n*(param.n0>0);
param.SO40 = 1.134*ones(param.N,1); %unit: 10^-10 mol
param.S20 = zeros(param.N,1);
param.Pb20 = 2.26*ones(param.N,1); %unit: 10^-10 mol
param.PbS0 = zeros(param.N,1);

