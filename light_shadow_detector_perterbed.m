% 20.320 Final Project
% Po-Han Lin
% MATLAB 2018b
% Dec 2018

% Initial Conditions for the Light Model
x0 = [33.9; 8510; 37819; 81.8; 366.8;...
      22.2; 4519.5; 8.9; 770.9; 128.8;...
      0.19; 5287.7; 0.0; 5.0; 0.4; 0];
x0_01 = [33.9; 8510; 37819; 81.8; 366.8;...
      22.2; 4519.5; 8.9; 770.9; 128.8;...
      0.19; 5287.7; 0.0; 5.0; 0.4; 0.1];
x0_1 = [33.9; 8510; 37819; 81.8; 366.8;...
      22.2; 4519.5; 8.9; 770.9; 128.8;...
      0.19; 5287.7; 0.0; 5.0; 0.4; 1];
x0_10 = [33.9; 8510; 37819; 81.8; 366.8;...
      22.2; 4519.5; 8.9; 770.9; 128.8;...
      0.19; 5287.7; 0.0; 5.0; 0.4; 10];
x0_20 = [33.9; 8510; 37819; 81.8; 366.8;...
      22.2; 4519.5; 8.9; 770.9; 128.8;...
      0.19; 5287.7; 0.0; 5.0; 0.4; 20];
x0_30 = [33.9; 8510; 37819; 81.8; 366.8;...
      22.2; 4519.5; 8.9; 770.9; 128.8;...
      0.19; 5287.7; 0.0; 5.0; 0.4; 30];
x0_40 = [33.9; 8510; 37819; 81.8; 366.8;...
      22.2; 4519.5; 8.9; 770.9; 128.8;...
      0.19; 5287.7; 0.0; 5.0; 0.4; 40];


  
%% Sensitivity curve  
  
b1 =  0.0426;
b2 =  143.6;
b3 =  0.129;

hv = 0:500;

exp = 1;
y1_K0  = b1 .* hv.^exp ./ (b2.^exp + hv.^exp);

exp = 2;
y2_K0  = b1 .* hv.^exp ./ (b2.^exp + hv.^exp);

exp = 3;
y3_K0  = b1 .* hv.^exp ./ (b2.^exp + hv.^exp);


figure 

plot(hv, y1_K0,'LineWidth',2)
hold on
plot(hv, y2_K0,'LineWidth',2)
hold on
plot(hv, y3_K0,'LineWidth',2)
legend('exp = 1', 'exp = 2', 'exp = 3')
ylabel('Photosynthesis rate (molecule/min)')
xlabel('Light intensity (hv)')
set(gca,'FontSize',12)


%% Plot glucose concentration vs. time  
tspan = 1:2000;

% Change exponent only
[t,x1_430] = ode15s(@(t,x) odefun(t,x,'16hr',430,1,0), tspan, x0);
[t,x1_100] = ode15s(@(t,x) odefun(t,x,'16hr',100,1,0), tspan, x0);
[t,x1_42]  = ode15s(@(t,x) odefun(t,x,'16hr',42,1,0), tspan, x0); 

[t,x2_430] = ode15s(@(t,x) odefun(t,x,'16hr',430), tspan, x0);
[t,x2_100] = ode15s(@(t,x) odefun(t,x,'16hr',100), tspan, x0);
[t,x2_42]  = ode15s(@(t,x) odefun(t,x,'16hr',42), tspan, x0);  

[t,x3_430] = ode15s(@(t,x) odefun(t,x,'16hr',430,3,0), tspan, x0);
[t,x3_100] = ode15s(@(t,x) odefun(t,x,'16hr',100,3,0), tspan, x0);
[t,x3_42]  = ode15s(@(t,x) odefun(t,x,'16hr',42,3,0), tspan, x0); 

t = 1:2000;

figure(1)
for i = 1:16
subplot(4,4,i)
plot(t,x1_42(:,i),'LineWidth',2); 
hold on
plot(t,x2_42(:,i),'LineWidth',2); 
hold on
plot(t,x3_42(:,i),'LineWidth',2);
%legend('1','2','3')
end

figure(2)
for i = 1:16
subplot(4,4,i)
plot(t,x1_100(:,i),'LineWidth',2); 
hold on
plot(t,x2_100(:,i),'LineWidth',2); 
hold on
plot(t,x3_100(:,i),'LineWidth',2);
%legend('1','2','3')
end

figure(3)
for i = 1:16
subplot(4,4,i)
plot(t,x1_430(:,i),'LineWidth',2); 
hold on
plot(t,x2_430(:,i),'LineWidth',2); 
hold on
plot(t,x3_430(:,i),'LineWidth',2);
%legend('1','2','3')
end


%% Plot peak values vs. light intensity

x1_neg_spikes1 = [];
x1_neg_spikes2 = [];
x1_neg_spikes3 = [];
x8_neg_spikes1 = [];
x8_neg_spikes2 = [];
x8_neg_spikes3 = [];
x1_pos_spikes1 = [];
x1_pos_spikes2 = [];
x1_pos_spikes3 = [];
x8_pos_spikes1 = [];
x8_pos_spikes2 = [];
x8_pos_spikes3 = [];

for ii = 0:200
    para = 2*ii;
    [t,x] = ode15s(@(t,x) odefun(t,x,'16hr',para,1,0), tspan, x0);
    x1_pos_spike = max(x(:,1));
    x1_neg_spike = min(x(:,1));
    x8_pos_spike = max(x(:,8));
    x8_neg_spike = min(x(:,8));
    x1_neg_spikes1 = [x1_neg_spikes1 x1_neg_spike];
    x8_neg_spikes1 = [x8_neg_spikes1 x8_neg_spike];
    x1_pos_spikes1 = [x1_pos_spikes1 x1_pos_spike];
    x8_pos_spikes1 = [x8_pos_spikes1 x8_pos_spike];
    
    
    [t,x] = ode15s(@(t,x) odefun(t,x,'16hr',para,2,0), tspan, x0);
    x1_pos_spike = max(x(:,1));
    x1_neg_spike = min(x(:,1));
    x8_pos_spike = max(x(:,8));
    x8_neg_spike = min(x(:,8));
    x1_neg_spikes2 = [x1_neg_spikes2 x1_neg_spike];
    x8_neg_spikes2 = [x8_neg_spikes2 x8_neg_spike];   
    x1_pos_spikes2 = [x1_pos_spikes2 x1_pos_spike];
    x8_pos_spikes2 = [x8_pos_spikes2 x8_pos_spike];
    
    [t,x] = ode15s(@(t,x) odefun(t,x,'16hr',para,3,0), tspan, x0);
    x1_pos_spike = max(x(:,1));
    x1_neg_spike = min(x(:,1));
    x8_pos_spike = max(x(:,8));
    x8_neg_spike = min(x(:,8));
    x1_neg_spikes3 = [x1_neg_spikes3 x1_neg_spike];
    x8_neg_spikes3 = [x8_neg_spikes3 x8_neg_spike];    
    x1_pos_spikes3 = [x1_pos_spikes3 x1_pos_spike];
    x8_pos_spikes3 = [x8_pos_spikes3 x8_pos_spike];
end



figure
subplot(1,2,1)
plot(0:2:400,x1_neg_spikes1,'LineWidth',2)
hold on
plot(0:2:400,x1_neg_spikes2,'LineWidth',2)
hold on
plot(0:2:400,x1_neg_spikes3,'LineWidth',2)
hold on
plot(0:2:400,x1_pos_spikes1,'LineWidth',2)
hold on
plot(0:2:400,x1_pos_spikes2,'LineWidth',2)
hold on
plot(0:2:400,x1_pos_spikes3,'LineWidth',2)
ylabel('Amount (molecule)')
xlabel('Light intensity hv')
title('free AtRGS1')
legend('neg1','neg2','neg3','pos1','pos2','pos3')
set(gca,'FontSize',12)

subplot(1,2,2)
plot(0:2:400,x8_neg_spikes1,'LineWidth',2)
hold on
plot(0:2:400,x8_neg_spikes2,'LineWidth',2)
hold on
plot(0:2:400,x8_neg_spikes3,'LineWidth',2)
hold on
plot(0:2:400,x8_pos_spikes1,'LineWidth',2)
hold on
plot(0:2:400,x8_pos_spikes2,'LineWidth',2)
hold on
plot(0:2:400,x8_pos_spikes3,'LineWidth',2)
ylabel('Amount (molecule)')
xlabel('Light intensity hv')
title('G\alpha^{GDP}')
legend('neg1','neg2','neg3','pos1','pos2','pos3')
set(gca,'FontSize',12)



%% Change ss value only
tspan = 1:8000;

[t,x_k0] = ode15s(@(t,x) odefun(t,x,'16hrlong',430,2,0), tspan, x0);
[t,x_k01] = ode15s(@(t,x) odefun(t,x,'16hrlong',430,2,0.1), tspan, x0_01);
[t,x_k1] = ode15s(@(t,x) odefun(t,x,'16hrlong',430,2,1), tspan, x0_1);
[t,x_k10] = ode15s(@(t,x) odefun(t,x,'16hrlong',430,2,10), tspan, x0_10);
[t,x_k20]  = ode15s(@(t,x) odefun(t,x,'16hrlong',430,2,20), tspan, x0_20);
[t,x_k30]  = ode15s(@(t,x) odefun(t,x,'16hrlong',430,2,30), tspan, x0_30);
[t,x_k40]  = ode15s(@(t,x) odefun(t,x,'16hrlong',430,2,40), tspan, x0_40);

t = 1:8000;

figure
for i = 1:16
subplot(4,4,i)
plot(t,x_k0(:,i),'LineWidth',2); 
hold on
plot(t,x_k01(:,i),'LineWidth',2); 
hold on
plot(t,x_k1(:,i),'LineWidth',2);
% hold on
% plot(t,x_k10(:,i),'LineWidth',2); 
% hold on
% plot(t,x_k20(:,i),'LineWidth',2); 
% hold on
% plot(t,x_k30(:,i),'LineWidth',2); 
% hold on
% plot(t,x_k40(:,i),'LineWidth',2);
% legend('x_k0','x_k10','x_k20','x_k30','x_k40')]
xlim([6000 8000])
xlabel('Time (min)')
ylabel('Amount (molecule)')
set(gca,'FontSize',12) 
end


figure
plot(t,100*percent_internal(x_k0),'LineWidth',2); 
hold on
plot(t,100*percent_internal(x_k01),'LineWidth',2); 
hold on
plot(t,100*percent_internal(x_k1),'LineWidth',2);
xlim([6000 8000])
ylim([0 100])
ylabel('% of internalized AtRGS1')
xlabel('Time (min)')
set(gca,'FontSize',12)




%% Change b3 value only
% [t,xx] = ode15s(@(t,x) odefun(t,x,'16hr',430,2,0,0.0129), tspan, x0);
% [t,xxx] = ode15s(@(t,x) odefun(t,x,'16hr',430,2,0,0.129), tspan, x0);
% [t,xxxx] = ode15s(@(t,x) odefun(t,x,'16hr',430,2,0,1.29), tspan, x0);
% 
% figure
% plot(t,xx(:,1),'LineWidth',1); 
% hold on
% plot(t,xxx(:,1),'LineWidth',1); 
% hold on
% plot(t,xxxx(:,1),'LineWidth',1); 


%% Light ODE Model

function dxdt = odefun(t,x,model,para,exp, K, new_b3)

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


if nargin == 5
    exp = exp;
    K = 0;
elseif nargin == 6
    exp = exp;
    K = K;
elseif nargin == 7
    exp = exp;
    K = K;
    b3 = new_b3;
else
    exp = 2;
    K = 0;
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


if strcmp(model,'16hrlong')
    if (6500 < t) && (t < 7500)
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


dxdt(16,:) = b1 * hv^exp / (b2^exp + hv^exp) - b3 * (x(16)-K);

end


% Calculating the percentage of internalized AtRGS1
function percentage = percent_internal(x)
percentage = x(:,12) ./ ( x(:,1) + x(:,2) + x(:,3) + x(:,4) + x(:,5) + x(:,11) + x(:,12) );
end


