function [datan, dataMet, dataMet5, dataS2, dataPbS, tdata]=spec_wrapper_function_dian(param)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters and Initial Distributions   %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parameters for simulation
N=param.N;

% initialize the distributions
n_max=param.n_max;
n0=param.n0;
Met0=param.Met0;
Met50=param.Met50;
SO40=param.SO40;
S20=param.S20;
Pb20=param.Pb20;
PbS0=param.PbS0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Preparation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tplot =1; %clf, drawnow, set(gcf,'renderer','zbuffer')
plotgap = round(tplot/param.dt);
dt = tplot/plotgap;
nplots = round(param.tmax/tplot);

t = 0;

options = odeset('InitialStep',0.1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% allocate matrices that store the dynamic fields %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

datan=zeros(nplots+1,N); % cells
datan(1,:)=n0;

dataMet=zeros(nplots+1,N); % Met
dataMet(1,:)=Met0;

dataMet5=zeros(nplots+1,N); % Met5
dataMet5(1,:)=Met50;

dataS2=zeros(nplots+1,N); % S2
dataS2(1,:)=S20;

dataPbS=zeros(nplots+1,N); % PbS
dataPbS(1,:)=PbS0;

tdata=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MAIN ITERATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:nplots
    for n=1:plotgap
        t=t+dt;
        
        inputs = [n0; Met0; Met50; S20; PbS0];
    
        sol=ode45(@adr_func_dian,[0 dt],inputs,options,param);
        vec=(deval(sol,dt));
    
        n0 = vec(1:N);
        Met0 = vec(N+1:2*N);
        Met50 = vec(2*N+1:3*N);
        S20 = vec(3*N+1:4*N);
        PbS0 = vec(4*N+1:5*N);
    end
    
    
    %% load into the history vector
    datan(i+1,:)=n0;
    dataMet(i+1,:)=Met0;
    dataMet5(i+1,:)=Met50;
    dataS2(i+1,:)=S20;
    dataPbS(i+1,:)=PbS0;

    tdata = [tdata t];
end

