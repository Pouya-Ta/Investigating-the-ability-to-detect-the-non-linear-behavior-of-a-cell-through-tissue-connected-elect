% Pouya Taghipour    -----> 9933014
% Amirhossein Mohebi -----> 9933057

clear all
close all
clc

%%The time-dependent behavior of the potassium conductance...
...can be simulated using a program that integrates the equation for 

f = 5; %Hz
t = 0:1/f * 1e+3; % ms
E_field = 1; % v/m
v = 10; % mv

%%This function(dXdT_n) can be integrated using MATLAB to simulate a voltage-clamp experiment. For example, we can simulate
...the response to a voltage-clamp experiment, with voltage set to v = 10 mV, and initial condition n = 0 with the
...following script:

v = 10; xo = 0; g_K = 10;
[t,x] = ode23s(@dXdT_n,[0 12],xo,[],v);
g = g_K*(x(:,1).^4);
figure()
plot(t,g)
xlabel('t (ms)');
ylabel('g_{K}');
set(gca,'fontsize',16)

%% Simulation of voltage clamp experiments for potassium current.
%% Data of Hodgkin and Huxley, J. Physiol, 116:500-544, 1952. Each %%
% data set corresponds to a different voltage clamp, as indicated.
% Time (in ms) is the first row, g_K is the second row.

% Data of v=10 mV: 
data2 = [0.1640    0.3330    0.5860    0.7540    1.1100    1.4900    2.0000    2.8400    4.1900    6.3800    8.8000   11.2000
         0.2500    0.2600    0.2600    0.3100    0.3100    0.3200    0.4000    0.5000    0.7100    1.0000    1.3000    1.6000];

v = 10.001; xo = 0.28; g_K = 10; % voltage, condition, and conductivity value
 
[t,x] = ode23s(@dXdT_n,[0 12],xo,[],v);
g = g_K*(x(:,1).^4);
figure()
plot(t,g,'k-','linewidth',1.5); hold on;
plot(data2(1,:),data2(2,:),'ko','linewidth',1.5,'Markerfacecolor',[1 1 1],'Markersize',5); hold off
set(gca,'Fontsize',9,'Ytick',[2],'Xtick',[]);
box off;
axis([0 12 0 2.4]);
text(12.2,0.45*2.4,'$v = 10$~mV','interpreter','latex','fontsize',10);


%%This function(dXdT_mh) can be integrated using MATLAB to simulate a voltage-clamp experiment. For example, we can simulate...
...the response to a voltage-clamp experiment, with voltage set to v = 10 mV, and initial condition n = 0 with the
...following script:

v = 10; xo = [0 1]; g_Na = 10;
[t,x] = ode23s(@dXdT_mh,[0 12],xo,[],v);
g = g_Na*(x(:,1).^3).*x(:,2);
figure()
plot(t,g)
xlabel('t (ms)');
ylabel('g_{Na}');
set(gca,'fontsize',16)

%% Simulation of voltage clamp experiments for sodium current.
% Data of Hodgkin and Huxley, J. Physiol, 116:500-544, 1952. Each
% data set corresponds to a different voltage clamp, as indicated.
% Time (in ms) is the first row, g_Na is the second row.

% Data of  v = 10 mV:  
data2 = [0.1490    0.3400    0.5100    0.7010    1.1000    1.4900    1.9700    2.8000    4.1600    6.4100    8.7700   11.3000
         0.0400    0.0900    0.1100    0.1200    0.1200    0.1200    0.1100    0.1100    0.1000    0.0900    0.0800    0.0800];
		 
v = 10.0; xo = [0.05 0.7];
 
[t,x] = ode15s(@dXdT_mh,[0 12],xo,[],v);
g = g_Na*(x(:,1).^3).*x(:,2);
figure()
plot(t,g,'k-','linewidth',1.5); hold on;
plot(data2(1,:),data2(2,:),'ko','linewidth',1.5,'Markerfacecolor',[1 1 1],'Markersize',5); hold off
set(gca,'Fontsize',9,'Ytick',[0.1],'Xtick',[]);
box off; axis([0 12 0 0.15]);
text(12.2,0.30*0.15,'$v = 10$~mV','interpreter','latex','fontsize',10);



%% Script below simulates the Hodgkin-Huxley model with applied
%  current to generate Figure 8.12.
% Generate initial condition xo for simulation:
% add zero applied current:
I_app = 0;
[~,x] = ode15s(@dXdT_HH,[0 30],[0 0 0 0],[],I_app);
xo = x(end,:);

% Now Adding nonzero applied current:
I_app = 6.2;
[t,x] = ode15s(@dXdT_HH,[0 30],xo,[],I_app);
% Plotting computed action potential

figure();
clf; 
hold on;
plot(t,x(:,1),'k-','linewidth',1.5);
set(gca,'Fontsize',18);
hold off;
box on; 
axis([0 30 -20 120]);
xlabel('$t$ (ms)','interpreter','latex','fontsize',20);
ylabel('$v$ (mV)','interpreter','latex','fontsize',20);
% Obtaining and plotting conductivity values

for i = 1:length(x)
  [f, G3(i,:)] = dXdT_HH(0,x(i,:),I_app);
end

figure();
clf;
hold on;
plot(t,G3(:,1),'k','linewidth',1.5);
plot(t,G3(:,2),'k--','linewidth',1.5);
set(gca,'Fontsize',18);
hold off; 
box on; 
axis([0 30 0 45]);
xlabel('$t$ (ms)','interpreter','latex','fontsize',20);
ylabel('Conductance (mS$\cdot$cm$^{-2}$)','interpreter','latex','fontsize',20);
text(16.4,21,'$g_{Na}$','interpreter','latex','fontsize',20);
text(19,8,'$g_{K}$','interpreter','latex','fontsize',20);

%%function that returns the computed time derivative dx/dt as a function of n(t) and v
function [f] = dXdT_n(~,n,v)
% FUNCTION dXdT_N
% Inputs: t - time (milliseconds)
% x - vector of state variables {n}
% V - applied voltage (mV)
%
% Outputs: f - vector of time derivatives
% {dn/dt}
% alphas and betas:
a_n = 0.01*(10-v)./(exp((10-v)/10)-1);
b_n = 0.125*exp(-v/80);
% Computing derivative dn/dt:
f(1,:) = a_n*(1-n) - b_n*n;
end
%The time-dependent behavior of the sodium conductance can be simulated using a program that integrates the
...equations for 

function [f] = dXdT_mh(~,x,v)
% FUNCTION dXdT_MH
% Inputs: t - time (milliseconds)
% x - vector of state variables {m,h}
% V - applied voltage (mV)
%
% Outputs: f - vector of time derivatives
% {dm/dt,dh/dt}
% state variables
m = x(1);
h = x(2);
% alphas and betas:
a_m = 0.1*(25-v)/(exp((25-v)/10)-1);
b_m = 4*exp(-v/18);
a_h = 0.07*exp(-v/20);
b_h = 1 ./ (exp((30-v)/10) + 1);
% Computing derivatives:
f(1,:) = a_m*(1-m) - b_m*m;
f(2,:) = a_h*(1-h) - b_h*h;
end
%% Combining the equations for the gating variables with the equation for membrane potential,
... we can write a function for the full combined model:

function [f,varargout] = dXdT_HH(~,x,I_app)
% FUNCTION dXdT_HH
% Inputs: t - time (milliseconds)
% x - vector of state variables {v,m,n,h}
% I_app - applied current (microA cm^{-2})
%
% Outputs: f - vector of time derivatives
% {dv/dt,dm/dt,dn/dt,dh/dt}
% Resting potentials, conductivities, and capacitance:
V_Na = 115;
V_K = -12;
V_L = 10.6;
g_Na = 120;
g_K = 36;
g_L = 0.3;
C_m = 1e-6;
% State Variables:
v = x(1);
m = x(2);
n = x(3);
h = x(4);
% alphas and betas:
a_m = 0.1*(25-v)/(exp((25-v)/10)-1);
b_m = 4*exp(-v/18);
a_h = 0.07*exp(-v/20);
b_h = 1 ./ (exp((30-v)/10) + 1);
a_n = 0.01*(10-v)./(exp((10-v)/10)-1);
b_n = 0.125*exp(-v/80);
% Computing currents:
I_Na = (m^3)*h*g_Na*(v-V_Na);
I_K = (n^4)*g_K*(v-V_K);
I_L = g_L*(v-V_L);
% Computing derivatives:
f(1) = (-I_Na - I_K - I_L + I_app)/C_m;
f(2,:) = a_m*(1-m) - b_m*m;
f(3) = a_n*(1-n) - b_n*n;
f(4) = a_h*(1-h) - b_h*h;
% Outputting the conductivities
varargout{1} = [(m^3)*h*g_Na (n^4)*g_K g_L];
end
