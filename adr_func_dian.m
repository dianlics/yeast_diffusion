function out=adr_func_dian(~,vec,param)
    % parameters for simulation
    N = param.N;
    length = param.L;
    h = length/N;

    % cell growth
    dn = param.dn;
    n_max = param.n_max;
    miu_max = param.miu_max;
    K_Met = param.K_Met;

    % Met
    dMet = param.dMet;
    YMet_n = param.YMet_n;
    vmax_Met = param.vmax_Met;

    % inhibition of cell growth from Pb2+
    KI_Pb = param.KI_Pb;
    mI_Pb = param.mI_Pb;

    % inhibition of cell growth from S2-
    KI_S2 = param.KI_S2;
    mI_S2 = param.mI_S2;

    % Met5 activity
    YMet5_n=param.YMet5_n;
    KI_Met=param.KI_Met;
    mI_Met=param.mI_Met;
    vmax_Met5=param.vmax_Met5;
    delta_Met5=param.delta_Met5;
    
    % S2
    dS2 = param.dS2;
    Ksp = param.Ksp;
    vmax_S2 = param.vmax_S2;

    n = vec(1:N);
    Met = vec(N+1:2*N);
    Met5 = vec(2*N+1:3*N);
    S2 = vec(3*N+1:4*N);
    Pb2 = param.Pb20;

    %% forward and backward (no-flux condition)
    np=[n(2:N); n(N-1)];
    nm=[n(2); n(1:N-1)];

    Metp=[Met(2:N); Met(N-1)];
    Metm=[Met(2); Met(1:N-1)];

    Met5p=[Met5(2:N); Met5(N-1)];
    Met5m=[Met5(2); Met5(1:N-1)];

    S2p=[S2(2:N); S2(N-1)];
    S2m=[S2(2); S2(1:N-1)];

    rad_vec=linspace(0,length,N)';

    phi=ones(N,1).*(n>1); % make sure cell number bigger than 1

    %% advection
    dMet5dt_a = dn./(n+1e-5).*1/(4*h^2).*(Met5p-Met5m).*(np-nm);

    %% diffusion
    dndt_d=dn/h^2*(np-2*n+nm)+ dn/2/h./rad_vec.*(np-nm);
    dndt_d(1)=2*dn/h^2*(np(1)-2*n(1)+nm(1)); % origin is special: no 1/r term (singular), but factor 2

    dMetdt_d=dMet/h^2*(Metp-2*Met+Metm)+ dMet/2/h./rad_vec.*(Metp-Metm);
    dMetdt_d(1)=2*dMet/h^2*(Metp(1)-2*Met(1)+Metm(1)); % origin is special: no 1/r term (singular), but factor 2

    dS2dt_d=dS2/h^2*(S2p-2*S2+S2m)+ dS2/2/h./rad_vec.*(S2p-S2m);
    dS2dt_d(1)=2*dS2/h^2*(S2p(1)-2*S2(1)+S2m(1)); % origin is special: no 1/r term (singular), but factor 2

    %% reaction
    % miu
    %miu = miu_max*Met./ (K_Met+Met);
    miu = miu_max*Met./ (K_Met+Met)./ (1+(S2/KI_S2).^mI_S2);
    %miu = miu_max*Met./ (K_Met+Met)./ (1+(S2/KI_S2).^mI_S2)./ (1+(Pb2/KI_Pb).^mI_Pb);

    % cell
    dndt_r = miu.*n.*(1-n/n_max);

    % Met
    dMetdt_r = vmax_Met*n.*phi - YMet_n*miu.*n.*(1-n/n_max);

    % Met5
    dMet5dt_r = -Met5.*miu.*(1-n/n_max) + vmax_Met5./(1+(Met/KI_Met).^mI_Met).*phi - delta_Met5*Met5;

    % S2
    dS2dt_r = vmax_S2*Met5.*n.*phi - Ksp*S2;

    %%% Integrating three parts together
    dndt = dndt_d+dndt_r;
    dMetdt = dMetdt_d+dMetdt_r;
    dMet5dt = dMet5dt_a+dMet5dt_r;
    dS2dt = dS2dt_d+dS2dt_r;

    out = [dndt; dMetdt; dMet5dt; dS2dt];

