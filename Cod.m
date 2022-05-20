clear all
close all
clf

%% ----- Parameters-----%%
 
param.n=40;                                 %[cells] - Number of grid cells/weight groups/weight classes

param.w_inf= 70000;                         %[g] - weight where they only reproduce (same as depth)
param.w_offspring=0.006;                    %[g] - weight of the smallest ones
param.w_mature=31;                          %step they mature - when they begin to reproduce. Grid no. 36>2000g
param.psi_mature=zeros(1,param.n);          %vector with zeros 
param.psi_mature(param.w_mature:end)=1;     %point of maturation. A vector with 0's when not mature and 1's when mature. 
param.E=0.3;                                %[unitless]Reproductive efficiency
param.Rmax=300000000;                       %[recruits] - Maximum recruitment
param.A=10;                                 %[g^0.25/year]  Growth constant
param.a=0.3;                                %[unitless]     Mortality factor% Number of grid cells/points/size groups
param.F=1.8;
param.w=logspace(log10(param.w_offspring),log10(param.w_inf),param.n); %[1/g] - Grid definition. 
param.dw = gradient(param.w);                                          %grid size in g

%% ----- Initital conditions ----- %% 
N0=zeros(1,param.n);                        %[individuals] - Initial condition. Vector with zeros: No fish in any weight classes. 
N0(1)=1000000;                              %[individuals] - Initial condition. Number of recruits in the first weight class for the first year. 
%% ----- Run ODE23 -----%
[t,y] = ode23(@Cod_function,[0:30], N0, [], param);
%% ----- Make figures -----%

figure(1)
    contourf(t,param.w,real(log10(y))',0:0.3:10)
    set(gca,'yscale','log','FontName','Times New Roman','FontSize',14)
    c = colorbar;
    c.Label.String = ('log10 concentration[#]')
    title('Size distribution of North Atlantic Cod over 30 years')
    ylabel('weight [g]')
    xlabel('Time [years]')
    
    
figure(2) %Viser at den når equilibrium ved ca. år 30. 
plot(log10(y(end,1:end)),param.w,'.-','Linewidth',2)
    hold on
    plot(log10(y(29,1:end)),param.w,'--.','Linewidth',2)
    hold on
    plot(log10(y(20,1:end)),param.w,'--.','Linewidth',2)
    hold on
    plot(log10(y(10,1:end)),param.w,'--.','Linewidth',2)
    hold on
    plot(log10(y(5,1:end)),param.w,'--.','Linewidth',2)
    hold on
    plot(log10(y(3,1:end)),param.w,'--.','Linewidth',2)
    hold on
    plot(log10(y(2,1:end)),param.w,'--.','Linewidth',2)
    set(gca,'FontName','Times New Roman','FontSize',14)
    legend('year 30','year 29','year 20','year 10','year 5','year 3','year 2')
    ylabel('weight[g]')
    xlabel('log10 concentration[#]')
    xlim([0 10])
    
figure(4) %plot of natural mortality
    mu_natural=param.a*param.A*param.w.^-0.25;
    plot(real(log10(param.w)),mu_natural, 'Linewidth',2)
    set(gca,'FontName','Times New Roman','FontSize',14)
    ylabel('mu natural')
    xlabel('weight [g]')
    xlim([-2.2 5])
    

    
%the below is just y(end,:) with 1 heaviside
Fish_w_2heaviside=[1222557544.21028,791790947.584071,512776352.853599,332107947.120920,215006209.969823,139308480.321822,90121474.2393467,58418639.7258882,37779352.0217884,24481107.6320504,15837479.9750706,10255188.9847749,6636166.53754561,4294922.62381483,2779031.41429935,1797985.63679961,1163067.70595264,752223.472399510,486407.075156236,314450.515570447,203231.057173730,131309.695126006,84811.1894142183,54756.5739532491,35336.1120736069,22791.2221711952,14690.7012852632,9462.23510517697,6089.18065680199,3914.33851197559,1304.89830729411,262.813465622833,46.6821313726193,11.5772285556719,7.35860402383618,4.65344057905346,2.92029001445022,1.80945264824636,1.09255375650464,0.742317533512193]
%the below is just y(end,:) with 2 heaviside
Fish_w_heaviside=[1222557544.22363,791788441.272703,512790840.755383,332068339.557772,215074411.898637,139225360.349092,90197860.1621395,58363569.5884831,37811348.3920822,24465833.4580181,15843557.7444059,10253150.7614387,6636747.59405950,4294780.85073724,2779061.18030069,1797980.23607294,1163068.55564390,752223.356196051,486407.089005707,314450.514135319,203231.057306344,131309.695117070,84811.1894158320,54756.5739537869,35336.1120739840,22791.2221714368,14690.7012854190,9462.23510527734,6089.18065686658,3914.33851201712,1884.79602340434,676.533845365247,223.475497132064,66.8659325833561,17.7494638519045,4.06280819919423,0.769575376553043,0.113042691367441,0.0114397572723680,0.000752566573309423]    

figure(6)
plot(log10(Fish_w_heaviside),param.w,'r-','Linewidth',2)
hold on
plot(log10(Fish_w_2heaviside),param.w,'b-','Linewidth',2)
xlim([0 10])
legend('Fishing with 1 heaviside','Fishing with 2 heaviside')
    

    Fish_w_heaviside_bio=Fish_w_heaviside.*param.dw; %to get in biomass (N*w)
    Fish_w_heaviside_bio_total=Fish_w_heaviside_bio(:,31:40); %biomass only of the fished columns
    Fish_w_heaviside_bio_total=sum(Fish_w_heaviside_bio_total) %all fished columns summed together should give total biomass LEFT 
    
    Fish_w_2heaviside_bio=Fish_w_2heaviside.*param.dw;
    Fish_w_2heaviside_bio_total=Fish_w_2heaviside_bio(:,31:34);
    Fish_w_2heaviside_bio_total=sum(Fish_w_2heaviside_bio_total)