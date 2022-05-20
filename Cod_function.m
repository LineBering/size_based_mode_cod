function dydt =Cod_funtion(t,y,param)  

N=y(1:param.n);

%% -----Advection(grwth)-----%%

for i=2:param.n 
    Ja_Growth(i)=(param.A*param.w(i-1)^(3/4))*(1-(param.w(i-1)/param.w_inf)^(1/4))*N(i-1); 
   
    
end

Rp=sum(param.psi_mature.*param.A*(param.w_inf^-0.25).*N.*param.dw);

%boundary flux, R
Ja_Growth(1)= (param.E*Rp)/(param.Rmax+param.E*Rp)*param.Rmax; 

%boundary flux when they are infinate large is 0
Ja_Growth(param.n+1)=0; 

 for i=1:param.n 
     dNdt(i)=-(Ja_Growth(i+1)-Ja_Growth(i))/param.dw(i); 
 end

%% ----- mu -----%%
mu_natural=param.a*param.A*param.w.^-0.25; %function of natural mortality
mu_fishing=heaviside(param.w -2000)*param.F;
mu_fishing=mu_fishing-heaviside(param.w -6000)*param.F;
%% ----Update derivaties---%% 

dNdt=dNdt-mu_natural.*N'-mu_fishing.*N'; % our reactive term 

dydt=dNdt';
end