clear all
close all
clf
%% Grid Sensitivity Analysis
%In this script we show the sensitivity of the grids. Since dw is on a log
%scale, the grid sizes changes with the amount of n. We therefor run with
%different n's. 

for  i=[20:2:36] 

    param.A=10; 
    param.a=0.1; 
    param.n=i;                  
    param.w_inf= 70000;           
    param.w_mature=31;            
    param.w_offspring=0.006;       
    param.w=logspace(log10(param.w_offspring),log10(param.w_inf),param.n); 
    param.dw = gradient(param.w);
    param.psi_mature=zeros(1,param.n);
    param.psi_mature(param.w_mature:end)=1;
    param.E=0.3;
    param.Rmax=300000000;
    param.K=param.A*param.w_inf^-0.25; 
    param.F=0%1.8%5;
    N0=zeros(1,param.n);
    N0(1)=1000000;
    [t,y] = ode23(@Cod_function,[0:30], N0, [], param);

    plot(log10(y(end,:)),param.w,'--.', 'Linewidth',1)
    hold on
    drawnow
end
legend('n=20','n=22','n=24','n=26','n=28','n=30','n=32','n=34','n=36')
set(gca,'FontName','Times New Roman','FontSize',14)
title('Sensitivity to n')
xlabel('log10 concentration [#]')
ylabel('Weight [g]')

xlim([-25 10])
