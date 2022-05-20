clear all
close all
clf
%% Fishing mortality
%In this script we show the biomass compared to different fishing mortalities 

for  i=[0:0.1:6] 
    param.F=i
    param.A=10; 
    param.a=0.3; 
    param.n=40;                  
    param.w_inf= 50000;           
    param.w_mature=33;            
    param.w_offspring=0.006;     
    param.w=logspace(log10(param.w_offspring),log10(param.w_inf),param.n); 
    param.dw = gradient(param.w);
    param.psi_mature=zeros(1,param.n);
    param.psi_mature(param.w_mature:end)=1;
    param.E=0.3;
    param.Rmax=300000000;
    

    N0=zeros(1,param.n);
    N0(1)=1000000;
    [t,y] = ode23(@Cod_function,[0:30], N0, [], param);
    
    %36:40 added together BIOMASS(N*W)
     yy=y(end,36:40).*param.dw(:,36:40);    
     yyy=sum(yy); %N*W--> total biomass of fished weight classes 
     L=y(end,:).*param.dw; %total biomass of everything in our basin
     LL=sum(L);
     
  %plot
     yyaxis left
     plot(param.F,LL,'b*','Linewidth',2)
     ylabel('Biomass (N*w)')
     set(gca,'FontName','Times New Roman','FontSize',14)

     hold on
     drawnow
     yyaxis right
     ylabel('Biomass (N*w)')
     plot(param.F,yyy,'r*','Linewidth',2)
   
  

end

    xlabel('Fishing mortality')
    legend('Biomass of total cod population','Biomass of target group of cod')
 