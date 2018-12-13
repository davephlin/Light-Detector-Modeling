% 20.320 Final Project
% Po-Han Lin
% MATLAB 2018b
% Dec 2018

% Initial Condition for the Light Model
x0 = [33.9; 8510; 37819; 81.8; 366.8;...
      22.2; 4519.5; 8.9; 770.9; 128.8;...
      0.19; 5287.7; 0.0; 5.0; 0.4; 0];
    
    
%% Figure 2a
tspan = 1:2000;

% Solving ODEs
[t,x430] = ode15s(@(t,x) odefun(t,x,'16hr',430), tspan, x0);
[t,x100] = ode15s(@(t,x) odefun(t,x,'16hr',100), tspan, x0);
[t,x42]  = ode15s(@(t,x) odefun(t,x,'16hr',42), tspan, x0);

% Plotting the concentrations with respect to time
figure()

subplot(3,2,1)
plot(t,x430(:,1),'LineWidth',2);
hold on
plot(t,x100(:,1),'LineWidth',2);
hold on
plot(t,x42(:,1),'LineWidth',2);
title('free AtRGS1')
xlabel('Time (min)')
ylabel('Amount (molecule)')
lgd = legend('430','100','42');
title(lgd,'Light Intensity')
xlim([0 2000])
set(gca,'FontSize',12)

subplot(3,2,2)
plot(t,x430(:,4),'LineWidth',2);
hold on
plot(t,x100(:,4),'LineWidth',2);
hold on
plot(t,x42(:,4),'LineWidth',2);
title('AtRGS1:G\alpha^{GDP}')
xlabel('Time (min)')
ylabel('Amount (molecule)')
xlim([0 2000])
set(gca,'FontSize',12) 

subplot(3,2,3)
plot(t,x430(:,8),'LineWidth',2);
hold on
plot(t,x100(:,8),'LineWidth',2);
hold on
plot(t,x42(:,8),'LineWidth',2);
title('G\alpha^{GDP}')
xlabel('Time (min)')
ylabel('Amount (molecule)')
xlim([0 2000])
set(gca,'FontSize',12) 

subplot(3,2,4)
plot(t,x430(:,5),'LineWidth',2);
hold on
plot(t,x100(:,5),'LineWidth',2);
hold on
plot(t,x42(:,5),'LineWidth',2);
title('AtRGS1:G\alpha^{GTP}')
xlabel('Time (min)')
ylabel('Amount (molecule)')
xlim([0 2000])
set(gca,'FontSize',12) 

subplot(3,2,5)
plot(t,x430(:,9),'LineWidth',2);
hold on
plot(t,x100(:,9),'LineWidth',2);
hold on
plot(t,x42(:,9),'LineWidth',2);
title('G\alpha^{GTP}')
xlabel('Time (min)')
ylabel('Amount (molecule)')
xlim([0 2000])
set(gca,'FontSize',12) 

subplot(3,2,6)
plot(t,x430(:,10),'LineWidth',2);
hold on
plot(t,x100(:,10),'LineWidth',2);
hold on
plot(t,x42(:,10),'LineWidth',2);
title('G\beta\gamma')
xlabel('Time (min)')
ylabel('Amount (molecule)')
xlim([0 2000])
set(gca,'FontSize',12) 


%% Figure 2bc

tspan = 1:2000;

% Keep track of different key values
% pos_spike - the spike from dark to light at the 8th hour
% neg_spike - the spike from light to dark at the 16th hour
% ss - steady state value under light
x1_pos_spikes =[];
x1_neg_spikes = [];
x1_sss = [];
x8_pos_spikes =[];
x8_neg_spikes = [];
x8_sss = [];

% Solving ODEs with various irradiance intensity from 0 to 400
for ii = 0:200
    para = 2*ii;
    [t,x] = ode15s(@(t,x) odefun(t,x,'16hr',para), tspan, x0);
    x1_pos_spike = max(x(:,1));
    x1_neg_spike = min(x(:,1));
    x1_ss = x(1000,1);
    x8_pos_spike = max(x(:,8));
    x8_neg_spike = min(x(:,8));
    x8_ss = x(1000,8);
    
    x1_pos_spikes =[x1_pos_spikes x1_pos_spike];
    x1_neg_spikes = [x1_neg_spikes x1_neg_spike];
    x1_sss = [x1_sss x1_ss];
    x8_pos_spikes =[x8_pos_spikes x8_pos_spike];
    x8_neg_spikes = [x8_neg_spikes x8_neg_spike];
    x8_sss = [x8_sss x8_ss];
end

% Plotting the concentrations with respect to irradiance intensity
figure()

subplot(1,2,1)
plot(0:2:400,x1_pos_spikes,'LineWidth',2);
hold on
plot(0:2:400,x1_neg_spikes,'LineWidth',2);
hold on
plot(0:2:400,x1_sss,'LineWidth',2);
title('free AtRGS1')
xlabel('Irradiane Intensity')
ylabel('Amount (molecule)')
legend('pos','neg','ss')
set(gca,'FontSize',12) 
title('free AtRGS1')

subplot(1,2,2)
plot(0:2:400,x8_pos_spikes,'LineWidth',2);
hold on
plot(0:2:400,x8_neg_spikes,'LineWidth',2);
hold on
plot(0:2:400,x8_sss,'LineWidth',2);
title('free AtRGS1')
xlabel('Irradiane Intensity')
ylabel('Amount (molecule)')
legend('pos','neg','ss')
set(gca,'FontSize',12) 
title('G\alpha^{GDP}')

%% Figure 3a
tspan = 1:3500;

% Solving ODEs
[t,x430_2] = ode45(@(t,x) odefun(t,x,'sinfun',430), tspan, x0);
[t,x100_2] = ode45(@(t,x) odefun(t,x,'sinfun',100), tspan, x0);
[t,x42_2]  = ode45(@(t,x) odefun(t,x,'sinfun',42), tspan, x0);

% Plotting the concentrations with respect to time
figure

subplot(2,3,1)
plot(t,x430_2(:,1),'LineWidth',2);
hold on
plot(t,x100_2(:,1),'LineWidth',2);
hold on
plot(t,x42_2(:,1),'LineWidth',2);
title('free AtRGS1')
xlabel('Time (min)')
ylabel('Amount (molecule)')
lgd = legend('430','100','42');
title(lgd,'Light Intensity')
xlim([1750 3250])
ylim([0 50])
set(gca,'FontSize',12) 

subplot(2,3,2)
plot(t,x430_2(:,4),'LineWidth',2);
hold on
plot(t,x100_2(:,4),'LineWidth',2);
hold on
plot(t,x42_2(:,4),'LineWidth',2);
title('AtRGS1:G\alpha^{GDP}')
xlabel('Time (min)')
ylabel('Amount (molecule)')
xlim([1750 3250])
ylim([0 100])
set(gca,'FontSize',12) 

subplot(2,3,3)
plot(t,x430_2(:,5),'LineWidth',2);
hold on
plot(t,x100_2(:,5),'LineWidth',2);
hold on
plot(t,x42_2(:,5),'LineWidth',2);
title('AtRGS1:G\alpha^{GTP}')
xlabel('Time (min)')
ylabel('Amount (molecule)')
xlim([1750 3250])
set(gca,'FontSize',12) 

subplot(2,3,4)
plot(t,x430_2(:,8),'LineWidth',2);
hold on
plot(t,x100_2(:,8),'LineWidth',2);
hold on
plot(t,x42_2(:,8),'LineWidth',2);
title('G\alpha^{GDP}')
xlabel('Time (min)')
ylabel('Amount (molecule)')
xlim([1750 3250])
ylim([0 13])
set(gca,'FontSize',12) 


subplot(2,3,5)
plot(t,x430_2(:,9),'LineWidth',2);
hold on
plot(t,x100_2(:,9),'LineWidth',2);
hold on
plot(t,x42_2(:,9),'LineWidth',2);
title('G\alpha^{GTP}')
xlabel('Time (min)')
ylabel('Amount (molecule)')
xlim([1750 3250])
set(gca,'FontSize',12) 

subplot(2,3,6)
plot(t,x430_2(:,10),'LineWidth',2);
hold on
plot(t,x100_2(:,10),'LineWidth',2);
hold on
plot(t,x42_2(:,10),'LineWidth',2);
title('G\beta\gamma')
xlabel('Time (min)')
ylabel('Amount (molecule)')
xlim([1750 3250])
set(gca,'FontSize',12) 

%% Figure 3bcd

tspan = 1:3500;

% Solving ODEs
[t,xtri10_0] = ode45(@(t,x) odefun(t,x,'tri10',0), tspan, x0);
[t,xtri10_5] = ode45(@(t,x) odefun(t,x,'tri10',5), tspan, x0);
[t,xtri10_10] = ode45(@(t,x) odefun(t,x,'tri10',10), tspan, x0);
[t,xtri10_15] = ode45(@(t,x) odefun(t,x,'tri10',15), tspan, x0);

[t,xtri60_0] = ode45(@(t,x) odefun(t,x,'tri60',0), tspan, x0);
[t,xtri60_5] = ode45(@(t,x) odefun(t,x,'tri60',5), tspan, x0);
[t,xtri60_10] = ode45(@(t,x) odefun(t,x,'tri60',10), tspan, x0);
[t,xtri60_15] = ode45(@(t,x) odefun(t,x,'tri60',15), tspan, x0);

[t,xtri60b_0] = ode45(@(t,x) odefun(t,x,'tri60base',0), tspan, x0);
[t,xtri60b_5] = ode45(@(t,x) odefun(t,x,'tri60base',5), tspan, x0);
[t,xtri60b_10] = ode45(@(t,x) odefun(t,x,'tri60base',10), tspan, x0);
[t,xtri60b_15] = ode45(@(t,x) odefun(t,x,'tri60base',15), tspan, x0);


% Plotting the concentrations with respect to time
figure

subplot(1,3,1)
plot(t,xtri10_0(:,1),'LineWidth',2);
hold on
plot(t,xtri10_5(:,1),'LineWidth',2);
hold on
plot(t,xtri10_10(:,1),'LineWidth',2);
hold on
plot(t,xtri10_15(:,1),'LineWidth',2);
title('free AtRGS1')
xlabel('Time (min)')
ylabel('Amount (molecule)')
xlim([1750 3250])
set(gca,'FontSize',12)
legend('0','5','10','15')


subplot(1,3,2)
plot(t,xtri60_0(:,1),'LineWidth',2);
hold on
plot(t,xtri60_5(:,1),'LineWidth',2);
hold on
plot(t,xtri60_10(:,1),'LineWidth',2);
hold on
plot(t,xtri60_15(:,1),'LineWidth',2);
title('free AtRGS1')
xlabel('Time (min)')
ylabel('Amount (molecule)')
xlim([1750 3250])
set(gca,'FontSize',12)
legend('0','5','10','15')


subplot(1,3,3)
plot(t,xtri60b_0(:,1),'LineWidth',2);
hold on
plot(t,xtri60b_5(:,1),'LineWidth',2);
hold on
plot(t,xtri60b_10(:,1),'LineWidth',2);
hold on
plot(t,xtri60b_15(:,1),'LineWidth',2);
title('free AtRGS1')
xlabel('Time (min)')
ylabel('Amount (molecule)')
xlim([1750 3250])
set(gca,'FontSize',12)
legend('0','5','10','15')

%% Figure 4ab

tspan = 1:250;
[t,x2] = ode45(@(t,x) odefun(t,x,'2min',430), tspan, x0);
[t,x4] = ode45(@(t,x) odefun(t,x,'4min',430), tspan, x0);
[t,x10] = ode45(@(t,x) odefun(t,x,'10min',430), tspan, x0);
[t,xc] = ode45(@(t,x) odefun(t,x,'const',430), tspan, x0);

figure

subplot(1,2,1)
plot(t,x2(:,1),'LineWidth',2);
hold on
plot(t,x4(:,1),'LineWidth',2);
hold on
plot(t,x10(:,1),'LineWidth',2);
hold on
plot(t,xc(:,1),'LineWidth',2);
title('free AtRGS1')
xlabel('Time (min)')
ylabel('Amount (molecule)')
xlim([80 250])
set(gca,'FontSize',12)
legend('2 min', '4 min', '10 min', 'C')

subplot(1,2,2)
plot(t,100*percent_internal(x2),'LineWidth',2);
hold on
plot(t,100*percent_internal(x4),'LineWidth',2);
hold on
plot(t,100*percent_internal(x10),'LineWidth',2);
hold on
plot(t,100*percent_internal(xc),'LineWidth',2);
title('% of internalized AtRGS1')
xlabel('Time (min)')
ylabel('Amount (molecule)')
xlim([80 250])
ylim([0 100])
set(gca,'FontSize',12)
legend('2 min', '4 min', '10 min', 'C')

%% Figure 5

% Solving ODEs
tspan = 1:10000;
[t1,xwt] = ode45(@(t,x) odefun(t,x,'const500','WT'), tspan, x0);
[t1,xk1] = ode45(@(t,x) odefun(t,x,'const500','k1'), tspan, x0);
[t1,xk2] = ode45(@(t,x) odefun(t,x,'const500','k2'), tspan, x0);

tspan = 1:5000;
[t2,xxwt] = ode45(@(t,x) odefun(t,x,'sins','WT'), tspan, x0);
[t2,xxk1] = ode45(@(t,x) odefun(t,x,'sins','k1'), tspan, x0);
[t2,xxk2] = ode45(@(t,x) odefun(t,x,'sins','k2'), tspan, x0);
%%
% Plotting the concentrations with respect to time
figure

subplot(1,2,1)
plot(t1,percent_internal(xwt)*100,'LineWidth',2);
hold on
plot(t1,percent_internal(xk1)*100,'LineWidth',2);
hold on
plot(t1,percent_internal(xk2)*100,'LineWidth',2);

title('% of internalized AtRGS1')
xlabel('Time (min)')
ylabel('Amount (molecule)')
ylim([0 70])
set(gca,'FontSize',12)
legend('WT', '\Delta kinase 1', '\Delta kinase 2')

subplot(1,2,2)
plot(t2,percent_internal(xxwt)*100,'LineWidth',2);
hold on
plot(t2,percent_internal(xxk1)*100,'LineWidth',2);
hold on
plot(t2,percent_internal(xxk2)*100,'LineWidth',2);
title('% of internalized AtRGS1')
xlabel('Time (min)')
ylabel('Amount (molecule)')
ylim([0 50])
xlim([0 5000])
set(gca,'FontSize',12)
legend('WT', '\Delta kinase 1', '\Delta kinase 2')


%% Light ODE Model

function dxdt = odefun(t,x,model,para)

% rate constants
k1 =  6.21e-4;
k2 =  2.15;
k3 =  7.29e-3;
k4 =  1.89;
k5 =  4.95e-2;
k6 =  8.4;
k7 =  3.50e-2;
k8 =  4.91e-3;
k9 =  4.26e-6;
k10 = 4.65e-5;
k11 = 29.38;
k12 = 5.00;
k13 = 6.61e-4; 
k14 = 2;
k15 = 2.31e-1; 
k16 = 5.44e-3;
k17 = 5.09e-4;
k18 = 3.21e-2;
k19 = 1.52e-2;
k20 = 7.20e1;
k21 = 5.53e3;
k22 = 4.25e3;
k23 = 5.57e4;
k24 = 3.80e-3;
k25 = 3.04e-3;
k26 = 1.92;
k27 = 1.84e-1;
k28 = 3.48;
k29 = 4.66e4;
k30 = 2.41;
b1 =  0.0426;
b2 =  143.6;
b3 =  0.129;

dxdt = zeros(16,1);

% WT / dKinase1 / dKinase2 model
if strcmp(para,'k1')
    x(15) = 0;
elseif strcmp(para,'k2')
    x(14) = 0;
end

% ODE   
dxdt(1,:) =  - k17 * x(1) * x(7) - k18 * x(1) * x(9) - k12 * x(1) + k24 * x(5) + k25 * x(3) + k27 * x(12);
dxdt(2,:) =  k4 * x(3) + k8 * x(4) * x(10) - k6 * x(2) - k16 * x(2);
dxdt(3,:) =  k6 * x(2) - k4 * x(3) - k25 * x(3) + k17 * x(1) * x(7) + k13 * x(5) * x(10) - k11 * x(3) * (x(13)^k14 / (k26^k14 + x(13)^k14)) ;
dxdt(4,:) =  k4 * x(5) - k6 * x(4) - k8 * x(4) * x(10) + k16 * x(2);
dxdt(5,:) =  k6 * x(4) - k4 * x(5) - k3 * x(5) * (x(14) + x(15)) + k30 * x(11) + k18 * x(1) * x(9) - k2 * x(5) -k24 * x(5) - k13 * x(5) * x(10) + k11 * x(3) * (x(13)^k14 / (k26^k14 + x(13)^k14));
dxdt(6,:) =  k5 * x(7) + k7 * x(8) * x(10) - k6 * x(6) - k28 * x(6);
dxdt(7,:) =  k6 * x(6) - k5 * x(7) - k10 * x(7) - k17 * x(1) * x(7) + k9 * x(9) * x(10) + k25 * x(3);
dxdt(8,:) =  k5 * x(9) + k28 * x(6) - k6 * x(8) - k7 * x(8) * x(10);
dxdt(9,:) =  k20 * x(11) + k6 * x(8) + k10 * x(7) - k5 * x(9) - k18 * x(1) * x(9) - k9 * x(9) * x(10) + k2 * x(5) + k24 * x(5);
dxdt(10,:) = k10 * x(7) + k11 * x(3) * (x(13)^k14 / (k26^k14 + x(13)^k14)) - k8 * x(4) * x(10) - k7 * x(8) * x(10) - k9 * x(9) * x(10) + k28 * x(6) - k13 * x(5) * x(10) + k16 * x(2);
dxdt(11,:) = k3 * x(5) * (x(14) + x(15)) - k20 * x(11) - k30 * x(11);
dxdt(12,:) = k2 * x(5) + k12 * x(1) + k20 * x(11) - k27 * x(12);
dxdt(13,:) = k15 * (x(16) - x(13));
dxdt(14,:) = k1 * (k21 * x(10)^2 / (k22^2 + x(10)^2) - x(14));
dxdt(15,:) = k19 * (k23 * x(10)^2 / (k29^2 + x(10)^2) - x(15));

% Various light patterns
% Figure 2
if strcmp(model,'16hr')
    if (500 < t) && (t < 1500)
        hv = para;
    else
        hv = 0;
    end
end
% Figure 3a
if strcmp(model,'sinfun')
   if  (500 < t) && (t < 1500)
       hv = para*cos(abs(t-1000)/2000*2*pi);
   elseif (2000 < t) && (t < 3000)
       hv = para*cos(abs(t-2500)/2000*2*pi);
   else
       hv = 0;
   end
end
% Figure 3b
if strcmp(model,'tri10')
   if  (500 < t) && (t < 1500)
       hv = para*floor((500-abs(t-1000))/10);
   elseif (2000 < t) && (t < 3000)
       hv = para*floor((500-abs(t-2500))/10);
   else
       hv = 0;
   end
end
% Figure 3c
if strcmp(model,'tri60')
   if  (500 < t) && (t < 1500)
       hv = para*floor((500-abs(t-1000))/60);
   elseif (2000 < t) && (t < 3000)
       hv = para*floor((500-abs(t-2500))/60);
   else
       hv = 0;
   end
end
% Figure 3d
if strcmp(model,'tri60base')
   if  (500 < t) && (t < 1500)
       hv = 50 + para*floor((500-abs(t-1000))/60);
   elseif (2000 < t) && (t < 3000)
       hv = 50 + para*floor((500-abs(t-2500))/60);
   else
       hv = 50;
   end
end
% Figure 4
if strcmp(model,'2min')
   if  (100 < t) && (t <= 102)
       hv = para*floor((1000-abs(t-1000))/60);
   else
       hv = 0;
   end
end

if strcmp(model,'4min')
   if  (100 < t) && (t <= 104)
       hv = para;
   else
       hv = 0;
   end
end

if strcmp(model,'10min')
   if  (100 < t) && (t <= 110)
       hv = para;
   else
       hv = 0;
   end
end

if strcmp(model,'const')
   if  100 < t
       hv = para;
   else
       hv = 0;
   end
end
% Figure 5a
if strcmp(model,'const500')
   if t > 500
       hv = 100;
   else
       hv = 0;
   end
end
% Figure 5b
if strcmp(model,'sins')
   if (500 < t) && (t < 1500)
       hv = 100*cos(abs(t-1000)/2000*2*pi);
   elseif (2000 < t) && (t < 3000)
       hv = 100*cos(abs(t-2500)/2000*2*pi);
   elseif (3500 < t) && (t < 4500)
       hv = 100*cos(abs(t-4000)/2000*2*pi);
   else
       hv = 0;
   end
end


dxdt(16,:) = b1 * hv^2 / (b2^2 + hv^2) - b3 * x(16);

end


% Calculating the percentage of internalized AtRGS1
function percentage = percent_internal(x)
percentage = x(:,12) ./ ( x(:,1) + x(:,2) + x(:,3) + x(:,4) + x(:,5) + x(:,11) + x(:,12) );
end


