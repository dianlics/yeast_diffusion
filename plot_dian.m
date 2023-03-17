%% plotting data
datapos = linspace(0,param.L,param.N);

tplot = 1; %clf, drawnow, set(gcf,'renderer','zbuffer')
plotgap = round(tplot/param.dt);
dt = tplot/plotgap;
nplots = round(param.tmax/tplot);

t = 0;

for i=1:nplots
    figure(1)
    tcl = tiledlayout(3,2);
    
    nexttile
    plot(datapos,datan(i,:))
    xlim([0 param.L])
    ylim([0 param.n_max])
    %set(gca, 'YScale', 'log')
    title('cell denity')
    xlabel('radius (mm)')
    ylabel('cell density')
    
    nexttile
    plot(datapos,dataMet(i,:))
    xlim([0 param.L])
    %set(gca, 'YScale', 'log')
    title('Met concentration')
    xlabel('radius (mm)')
    ylabel('Met concentration')
    
    nexttile
    plot(datapos,dataMet5(i,:))
    xlim([0 param.L])
    %set(gca, 'YScale', 'log')
    title('Met5 proteins per cell concentration')
    xlabel('radius (mm)')
    ylabel('Met5 proteins per cell concentration')
    
    nexttile
    plot(datapos,dataMet5(i,:).*datan(i,:))
    xlim([0 param.L])
    %set(gca, 'YScale', 'log')
    title('total Met5 proteins concentration')
    xlabel('radius (mm)')
    ylabel('total Met5 proteins concentration')
    
    nexttile
    plot(datapos,dataS2(i,:))
    xlim([0 param.L])
    %set(gca, 'YScale', 'log')
    title('S2 concentration')
    xlabel('radius (mm)')
    ylabel('S2 concentration')

    nexttile
    plot(datapos,dataPbS(i,:).*datan(i,:))
    xlim([0 param.L])
    %set(gca, 'YScale', 'log')
    title('total PbS concentration')
    xlabel('radius (mm)')
    ylabel('total PbS concentration')

    title(tcl,['At time ',num2str(i),' h'])
end
